import logging
from typing import Optional


def setup_logger(log: Optional[logging.Logger] = None, log_path: Optional[str] = None, log_name: Optional[str] = None) -> logging.Logger:
    if log is None:
        log = logging.getLogger(log_name)
        log.setLevel(logging.DEBUG)
        log.propagate = False  # Prevent messages from being logged twice

    # Only set this once to avoid re-adding handlers on repeated calls
    if not any(isinstance(h, logging.StreamHandler) for h in log.handlers):
        stream_handler = logging.StreamHandler()
        stream_handler.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        stream_handler.setFormatter(formatter)
        log.addHandler(stream_handler)

    if log_path and not any(isinstance(h, logging.FileHandler) for h in log.handlers):
        file_handler = logging.FileHandler(log_path, mode='w')
        file_handler.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        log.addHandler(file_handler)

    return log
