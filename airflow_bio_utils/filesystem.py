
from __future__ import annotations
from fs import open_fs
import tempfile
from io import StringIO, BytesIO
import os
import uuid
from enum import Enum, unique
from typing import List, Optional
from pathlib import Path
from dataclasses import dataclass


def url_join_path(remote_url: str, *args: List[str]) -> str:
    if remote_url.count("//") == 0:
        return f"file:////{os.path.join(remote_url, *args)}"
    elif remote_url.count("//") == 1:
        remote_url = f"{remote_url}//{os.path.join(remote_url, *args)}"
    return os.path.join(remote_url, *args)


@unique
class ReadMode(str, Enum):
    READ = "r"
    READ_MODIFY = "r+"
    WRITE = "w"
    READ_BYTES = "rb"
    READ_MODIFY_BYTES = "rb+"
    WRITE_BYTES = "wb"


@dataclass
class FileURL:
    protocol: str
    filesystem_config: str
    path: str
    descriptor = None

    def __str__(self) -> str:
        return f"{self.protocol}://{self.filesystem_config}//{self.path}"

    @property
    def basename(self) -> str:
        if len(self.path) > 0:
            return os.path.basename(self.path)
        return f"{uuid.uuid4().hex}.tmp"

    @property
    def fs_url(self) -> str:
        return f"{self.protocol}://{self.filesystem_config}"

    @classmethod
    def parse_url(cls, remote_url: str) -> FileURL:
        if remote_url.count("//") == 0:
            remote_url = str(Path(remote_url).resolve())
            return FileURL(
                protocol="file",
                filesystem_config="",
                path=remote_url,
            )
        elif remote_url.count("//") == 1:
            print(remote_url)
            [protocol, filesystem_config] = remote_url.split("://")
            return FileURL(
                protocol=protocol,
                filesystem_config=filesystem_config,
                path="",
            )
        else:
            fs_url, _, fs_dir = remote_url.rpartition('//')
            [protocol, filesystem_config] = fs_url.split("://")
            return FileURL(
                protocol=protocol,
                filesystem_config=filesystem_config,
                path=fs_dir,
            )


@dataclass
class FileHandler:
    url: FileURL
    mode: ReadMode
    descriptor = None
    was_pulled: bool = False

    def pull(self) -> str:
        if not self.was_pulled:
            temp_folder = tempfile.mkdtemp()
            temp_filepath = os.path.join(temp_folder, self.url.basename)
            with self as fs:
                with open(temp_filepath, "w") as local_file:
                    local_file.write(fs.read())
            self.url = FileURL(
                protocol="file",
                filesystem_config="",
                path=temp_filepath,
            )
            self.was_pulled = True
            return temp_filepath
        else:
            return self.url.path

    @property
    def local_file_path(self) -> str:
        self.pull()
        return self.url.path

    @property
    def local_file_descriptor(self):
        return open(self.local_file_path, self.mode)

    @classmethod
    def open_fs(cls, url: str):
        if url.startswith("s3://"):
            bucket_name = url.replace("s3://", "").split("?")[0]
            from fs_s3fs import S3FS
            # Fix for S3FS, because we don't need to validate directories strictly
            # It speeds up everything
            return S3FS(bucket_name, strict=False)
        return open_fs(url)

    def __enter__(self):
        if self.url.protocol == "file":
            # Ensure path exists
            os.makedirs(os.path.dirname(self.url.path), exist_ok=True)
            self.descriptor = open(self.url.path, self.mode)
        else:
            if self.mode == ReadMode.WRITE:
                self.descriptor = StringIO("")
            elif self.mode == ReadMode.WRITE_BYTES:
                self.descriptor = BytesIO(b"")
            else:
                with FileHandler.open_fs(self.url.fs_url) as fs:
                    self.descriptor = StringIO(fs.readtext(self.url.path))
        return self.descriptor

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.url.protocol == "file":
            self.descriptor.close()
        else:
            if self.mode == ReadMode.WRITE:
                with FileHandler.open_fs(self.url.fs_url) as fs:
                    self.descriptor.seek(0)
                    if self.url.protocol == "s3":
                        # writetext(...) has bug with strict=False for S3FS
                        # so we use writebytes
                        fs.writebytes(self.url.path, self.descriptor.read().encode())
                    else:
                        fs.writetext(self.url.path, self.descriptor.read())
            elif self.mode == ReadMode.WRITE_BYTES:
                with FileHandler.open_fs(self.url.fs_url) as fs:
                    self.descriptor.seek(0)
                    fs.writebytes(self.url.path, self.descriptor.read())
            self.descriptor.close()


def open_url(url: str, mode: ReadMode = ReadMode.READ) -> FileHandler:
    parsed_url = FileURL.parse_url(url)
    return FileHandler(url=parsed_url, mode=mode)
