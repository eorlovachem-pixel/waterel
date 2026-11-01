#!/usr/bin/env python3
"""Non-destructive upload of staged Cl/Br folders to the shared Google Drive.

Requires a Desktop OAuth client JSON at GOOGLE_OAUTH_CLIENT or
~/hydration-research/google_oauth_client.json. Never overwrites a remote file
whose SHA-256 already matches. Does not upload fabricated QC.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STAGE = ROOT / "drive_staging"
PARENT_ID = "14ltUb63LmmihKGgzMf-geoXxPUmNLSp2"
TOKEN = Path.home() / "hydration-research" / "google_drive_token.json"
CLIENT_CANDIDATES = [
    Path(os.environ["GOOGLE_OAUTH_CLIENT"]) if os.environ.get("GOOGLE_OAUTH_CLIENT") else None,
    Path.home() / "hydration-research" / "google_oauth_client.json",
    ROOT / "google_oauth_client.json",
]
SCOPES = ["https://www.googleapis.com/auth/drive.file"]


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def client_secret() -> Path:
    for p in CLIENT_CANDIDATES:
        if p and p.is_file():
            return p
    raise FileNotFoundError(
        "No Google OAuth client JSON. Place a Desktop client at "
        "~/hydration-research/google_oauth_client.json or set GOOGLE_OAUTH_CLIENT. "
        "Staging trees are already in drive_staging/."
    )


def main() -> None:
    try:
        from google.auth.transport.requests import Request
        from google.oauth2.credentials import Credentials
        from google_auth_oauthlib.flow import InstalledAppFlow
        from googleapiclient.discovery import build
        from googleapiclient.http import MediaFileUpload
    except ImportError as exc:
        raise SystemExit(
            "google-api-python-client / google-auth-oauthlib not installed. "
            "Local staging is complete; install the libraries to upload.\n"
            f"{exc}"
        ) from exc

    secret = client_secret()
    creds = None
    if TOKEN.is_file():
        creds = Credentials.from_authorized_user_file(str(TOKEN), SCOPES)
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
        else:
            flow = InstalledAppFlow.from_client_secrets_file(str(secret), SCOPES)
            creds = flow.run_local_server(port=0)
        TOKEN.parent.mkdir(parents=True, exist_ok=True)
        TOKEN.write_text(creds.to_json(), encoding="utf-8")

    service = build("drive", "v3", credentials=creds)

    def find_child(parent_id: str, name: str):
        q = (
            f"'{parent_id}' in parents and name = '{name}' "
            "and trashed = false"
        )
        res = service.files().list(q=q, fields="files(id,name,md5Checksum,sha256Checksum)").execute()
        files = res.get("files") or []
        return files[0] if files else None

    def ensure_folder(parent_id: str, name: str) -> str:
        existing = find_child(parent_id, name)
        if existing:
            return existing["id"]
        meta = {
            "name": name,
            "mimeType": "application/vnd.google-apps.folder",
            "parents": [parent_id],
        }
        created = service.files().create(body=meta, fields="id").execute()
        return created["id"]

    report = {"parent": PARENT_ID, "folders": {}}
    for folder in ("Cl-_eorlovachem-pixel", "Br-_eorlovachem-pixel"):
        local = STAGE / folder
        if not local.is_dir():
            raise FileNotFoundError(local)
        fid = ensure_folder(PARENT_ID, folder)
        uploaded = []
        for path in sorted(p for p in local.rglob("*") if p.is_file()):
            rel = path.relative_to(local)
            parent = fid
            for part in rel.parts[:-1]:
                parent = ensure_folder(parent, part)
            name = rel.name
            remote = find_child(parent, name)
            digest = sha256(path)
            if remote:
                uploaded.append({
                    "path": str(rel.as_posix()),
                    "action": "exists_left_in_place",
                    "id": remote["id"],
                    "local_sha256": digest,
                })
                continue
            media = MediaFileUpload(str(path), resumable=False)
            created = service.files().create(
                body={"name": name, "parents": [parent]},
                media_body=media,
                fields="id",
            ).execute()
            uploaded.append({
                "path": str(rel.as_posix()),
                "action": "created",
                "id": created["id"],
                "local_sha256": digest,
            })
        report["folders"][folder] = {"id": fid, "files": uploaded}
        url = f"https://drive.google.com/drive/folders/{fid}"
        report["folders"][folder]["url"] = url
        (local / "DRIVE_UPLOAD_REPORT.json").write_text(
            json.dumps(report["folders"][folder], indent=2) + "\n", encoding="utf-8"
        )
        print(folder, url)
    (STAGE / "DRIVE_UPLOAD_REPORT.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
