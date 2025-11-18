import unittest
from unittest.mock import patch, Mock
from pipeline_util.fetch_mnemonic import fetch_reference_proteome_upid


class TestFetchReferenceProteome(unittest.TestCase):
    @patch("pipeline_util.fetch_mnemonic.requests.get")
    def test_successful_response(self, mock_get):
        # Mock JSON response
        mock_get.return_value = Mock(
            status_code=200,
            json=lambda: {
                "results": [
                    {
                        "proteomeType": "Reference and representative proteome",
                        "superkingdom": "Eukaryota",
                        "taxonomy": {"mnemonic": "HUMAN"},
                    }
                ]
            },
        )

        mnemonic = fetch_reference_proteome_upid("Homo sapiens")
        self.assertEqual(mnemonic, "HUMAN")

    @patch("pipeline_util.fetch_mnemonic.requests.get")
    def test_no_valid_proteome(self, mock_get):
        mock_get.return_value = Mock(status_code=200, json=lambda: {"results": []})
        mnemonic = fetch_reference_proteome_upid("Unknown organism")
        self.assertIsNone(mnemonic)

    @patch("pipeline_util.fetch_mnemonic.requests.get")
    def test_request_exception(self, mock_get):
        import requests
        mock_get.side_effect = requests.RequestException("API Error")
        mnemonic = fetch_reference_proteome_upid("Homo sapiens")
        self.assertIsNone(mnemonic)


if __name__ == "__main__":
    unittest.main()
