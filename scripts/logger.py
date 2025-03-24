import logging
import os
from enum import IntEnum
from datetime import datetime
from functools import wraps
from typing import Optional
from rich.progress import Progress, BarColumn, TimeElapsedColumn, TimeRemainingColumn, TextColumn
import traceback


class Logger:
    LOG_DIR = "logs"

    class Level(IntEnum):
        NOTSET = logging.NOTSET
        DEBUG = logging.DEBUG
        INFO = logging.INFO
        WARNING = logging.WARNING
        ERROR = logging.ERROR
        CRITICAL = logging.CRITICAL

    def __init__(self,
                 name: Optional[str] = None,
                 level: 'Logger.Level' = Level.INFO,
                 log_to_console: bool = True,
                 log_to_file: bool = True,
                 draw_progress: bool = False):

        self.logger = logging.getLogger(name)
        self.logger.setLevel(level)

        os.makedirs(self.LOG_DIR, exist_ok=True)
        log_filename = os.path.join(
            self.LOG_DIR, f"log_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.log"
        )

        formatter = logging.Formatter('[%(asctime)s - %(levelname)s]: %(message)s',
                                      '%Y-%m-%d %H:%M:%S')

        if not self.logger.handlers:
            if log_to_console:
                console_handler = logging.StreamHandler()
                console_handler.setFormatter(formatter)
                self.logger.addHandler(console_handler)

            if log_to_file:
                file_handler = logging.FileHandler(log_filename, mode='a', encoding='utf-8')
                file_handler.setFormatter(formatter)
                self.logger.addHandler(file_handler)

        if draw_progress:
            self.progress = Progress(
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TimeElapsedColumn(),
                TimeRemainingColumn()
            )
            self.progress.start()
        else:
            self.progress = None

    def get_logger(self):
        return self.logger

    def trace_call(self, func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            func_name = func.__qualname__
            self.logger.info(f"[START] {func_name}")
            try:
                result = func(*args, **kwargs)
                self.logger.info(f"[END] {func_name}")
                return result
            except Exception as e:
                self.logger.error(f"[ERROR in {func_name}]: {e}")
                self.logger.debug(traceback.format_exc())
                raise
        return wrapper

    def log_error(self, message: str, exception: Optional[Exception] = None):
        if exception:
            self.logger.error(f"{message} | Exception: {str(exception)}")
            self.logger.debug(traceback.format_exc())
        else:
            self.logger.error(message)

    def progress_task(self, description, total=100):
        if self.progress:
            task = self.progress.add_task(description, total=total)
            return task
        return None

    def advance_progress(self, task, advance=1):
        if self.progress and task is not None:
            self.progress.update(task, advance=advance)

    def finish_progress(self, task):
        if self.progress and task is not None:
            self.progress.update(task, completed=100)
            self.progress.stop()