import os
from dataclasses import dataclass


@dataclass(frozen=True)
class AppConfig:
    app_name: str = "LeishDomains"
    default_session_prefix: str = "SESSION"

    @property
    def cwd(self) -> str:
        return os.getcwd()


