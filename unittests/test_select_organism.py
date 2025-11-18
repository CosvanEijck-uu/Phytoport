import io
import unittest
from unittest.mock import patch, MagicMock, mock_open
import requests

from pipeline_util.select_organism import (
    fetch_reference_proteome_upid,
    download_gz_fasta,
    download_reference_proteome_for_organism,
    download_reference_proteome_for_upids,
)


class TestFetchReferenceProteomeUpid(unittest.TestCase):
    @patch("pipeline_util.select_organism.requests.get")
    def test_fetch_valid_proteome(self, mock_get):
        mock_response = MagicMock()
        mock_response.json.return_value = {
            "results": [{"id": "UP000005640", "superkingdom": ["Eukaryota"]}]
        }
        mock_response.raise_for_status = MagicMock()
        mock_get.return_value = mock_response

        upid = fetch_reference_proteome_upid("Homo sapiens")
        self.assertEqual(upid, "UP000005640")

    @patch("pipeline_util.select_organism.requests.get")
    def test_fetch_returns_none_for_viruses(self, mock_get):
        mock_response = MagicMock()
        mock_response.json.return_value = {
            "results": [{"id": "UP000001", "superkingdom": ["viruses"]}]
        }
        mock_response.raise_for_status = MagicMock()
        mock_get.return_value = mock_response

        upid = fetch_reference_proteome_upid("Some virus")
        self.assertIsNone(upid)

    @patch("pipeline_util.select_organism.requests.get")
    def test_fetch_handles_request_exception(self, mock_get):
        mock_get.side_effect = requests.RequestException("Network error")
        upid = fetch_reference_proteome_upid("Arabidopsis thaliana")
        self.assertIsNone(upid)


class TestDownloadGzFasta(unittest.TestCase):
    @patch("pipeline_util.select_organism.requests.get")
    @patch("pipeline_util.select_organism.open", new_callable=mock_open)
    @patch("pipeline_util.select_organism.gzip.open")
    @patch("pipeline_util.select_organism.shutil.copyfileobj")
    @patch("pipeline_util.select_organism.os.remove")
    def test_download_and_extract(
        self, mock_remove, mock_copyfile, mock_gzip_open, mock_open_file, mock_get
    ):
        mock_response = MagicMock()
        mock_response.iter_content = lambda chunk_size: [b"fake gz data"]
        mock_response.headers = {"content-length": "12"}
        mock_response.raise_for_status = MagicMock()
        mock_get.return_value.__enter__.return_value = mock_response

        mock_gzip_open.return_value.__enter__.return_value = io.BytesIO(
            b">protein1\nMSEQ..."
        )

        download_gz_fasta("UP000005640", "test_out.fasta")

        mock_open_file.assert_called()
        mock_copyfile.assert_called()
        mock_remove.assert_called_with("test_out.fasta.gz")

    @patch("pipeline_util.select_organism.os.path.exists", return_value=True)
    @patch("pipeline_util.select_organism.requests.get")
    def test_skips_if_file_exists(self, mock_get, mock_exists):
        download_gz_fasta("UP000005640", "test_out.fasta")
        mock_get.assert_not_called()


class TestDownloadReferenceProteomeForOrganism(unittest.TestCase):
    @patch("pipeline_util.select_organism.download_gz_fasta")
    @patch("pipeline_util.select_organism.fetch_reference_proteome_upid", return_value="UP000005640")
    def test_download_reference_proteome(self, mock_fetch, mock_download):
        download_reference_proteome_for_organism("Homo sapiens")
        mock_fetch.assert_called_once_with("Homo sapiens")
        mock_download.assert_called_once_with("UP000005640", "Homo_sapiens_ref_proteome.fasta")

    @patch("pipeline_util.select_organism.download_gz_fasta")
    @patch("pipeline_util.select_organism.fetch_reference_proteome_upid", return_value=None)
    def test_download_reference_proteome_none(self, mock_fetch, mock_download):
        download_reference_proteome_for_organism("Unknown species")
        mock_fetch.assert_called_once()
        mock_download.assert_not_called()


class TestDownloadReferenceProteomeForUpids(unittest.TestCase):
    @patch("pipeline_util.select_organism.download_gz_fasta")
    def test_download_multiple_upids_with_prefixes(self, mock_download):
        upids = ["UP000006548", "UP000005640"]
        prefixes = ["Arabidopsis", "Yeast"]

        download_reference_proteome_for_upids(upids, prefixes)

        mock_download.assert_any_call("UP000006548", "Arabidopsis_ref_proteome.fasta")
        mock_download.assert_any_call("UP000005640", "Yeast_ref_proteome.fasta")
        self.assertEqual(mock_download.call_count, 2)

    @patch("pipeline_util.select_organism.download_gz_fasta")
    def test_download_multiple_upids_with_fewer_prefixes(self, mock_download):
        upids = ["UP000006548", "UP000005640"]
        prefixes = ["Arabidopsis"]

        download_reference_proteome_for_upids(upids, prefixes)

        mock_download.assert_any_call("UP000006548", "Arabidopsis_ref_proteome.fasta")
        mock_download.assert_any_call("UP000005640", "UP000005640_ref_proteome.fasta")
        self.assertEqual(mock_download.call_count, 2)

    @patch("pipeline_util.select_organism.download_gz_fasta")
    def test_download_multiple_upids_without_prefixes(self, mock_download):
        upids = ["UP000006548", "UP000005640"]

        download_reference_proteome_for_upids(upids)

        mock_download.assert_any_call("UP000006548", "UP000006548_ref_proteome.fasta")
        mock_download.assert_any_call("UP000005640", "UP000005640_ref_proteome.fasta")
        self.assertEqual(mock_download.call_count, 2)


if __name__ == "__main__":
    unittest.main()
