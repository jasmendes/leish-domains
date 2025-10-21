import logging
import os
from typing import Optional


_LOGGER: Optional[logging.Logger] = None


def get_logger(name: str = "leishdomains") -> logging.Logger:
    global _LOGGER
    if _LOGGER is not None:
        return _LOGGER

    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)

    formatter = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    log_file = os.path.join(os.getcwd(), f"{name}.log")
    try:
        file_handler = logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    except Exception:
        # File logging is optional; ignore failures silently
        pass

    _LOGGER = logger
    return logger


