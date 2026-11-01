# Team Google Drive inventory

Parent: https://drive.google.com/drive/folders/14ltUb63LmmihKGgzMf-geoXxPUmNLSp2  
K GWH folder: https://drive.google.com/drive/folders/1SdWG8Lic11bysMid9j76cnchUhpCUNMs

**Not inventoried in this session.** No Desktop OAuth client JSON was present at `~/hydration-research/google_oauth_client.json`, and `google-api-python-client` is not installed. Team GWH files were **not** downloaded and **not** mixed into this kit.

Local staging (gitignored): `drive_staging/K+_eorlovachem-pixel/`, `Cl-_eorlovachem-pixel/`, `Br-_eorlovachem-pixel/`.

To upload after credentials exist:

```bash
pip install --user google-api-python-client google-auth-oauthlib
# place Desktop OAuth client JSON at ~/hydration-research/google_oauth_client.json
python3 workflow.py upload-drive
python3 workflow.py verify-delivery
```
