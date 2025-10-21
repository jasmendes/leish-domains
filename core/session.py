import os
from datetime import datetime
from typing import List

from .config import AppConfig
from .logger import get_logger


class SessionManager:
    def __init__(self, config: AppConfig | None = None):
        self.config = config or AppConfig()
        self.logger = get_logger()

    def save(self, file_paths: List[str], directory: str | None = None) -> str:
        if not file_paths:
            raise ValueError("No files to save in session")

        target_dir = directory or self.config.cwd
        os.makedirs(target_dir, exist_ok=True)

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        session_name = f"{self.config.default_session_prefix}_{timestamp}.txt"
        session_path = os.path.join(target_dir, session_name)

        with open(session_path, "w", encoding="utf-8") as f:
            for path in file_paths:
                f.write(f"{path}\n")

        self.logger.info("Session saved: %s", session_path)
        return session_path

    def load(self, session_file: str) -> List[str]:
        if not os.path.isfile(session_file):
            raise FileNotFoundError(session_file)

        loaded: List[str] = []
        with open(session_file, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line:
                    loaded.append(line)

        self.logger.info("Session loaded: %s (%d entries)", session_file, len(loaded))
        return loaded


