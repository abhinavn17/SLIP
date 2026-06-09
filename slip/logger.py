import sys
import os
import shutil
import glob
from datetime import datetime
import pytz
import atexit

class PrintLogger:
    """
    Redirects all print statements to a log file and optionally to terminal.
    When closed, moves all .log files from the current working directory
    into log_dir.
    """

    def __init__(self, name="slip_run", to_terminal=True, log_dir="logs"):
        """
        Initialize the logger.

        Parameters
        ----------
        name : str
            Base name for the log file.
        to_terminal : bool
            Whether to still print to terminal.
        log_dir : str
            Directory to store the log file.
        """
        self.to_terminal = to_terminal
        self.log_dir = log_dir
        self._closed = False

        # Create log directory if it doesn't exist
        os.makedirs(self.log_dir, exist_ok=True)

        # Get current time in Asia/Kolkata timezone
        tz_IN = pytz.timezone("Asia/Kolkata")
        timestamp = datetime.now(tz_IN).strftime("%Y-%m-%d_%H%M%S")

        self.log_file = os.path.join(
            self.log_dir,
            f"{name}_{timestamp}.log"
        )

        self.logfile = open(self.log_file, "a", encoding="utf-8")

        # Save original stdout
        self.terminal = sys.stdout if self.to_terminal else sys.__stdout__

        # Redirect stdout
        sys.stdout = self

        self.write(f"Logging started. Log file: {self.log_file}\n")
        atexit.register(self.close)

    def write(self, message):
        if self.terminal:
            self.terminal.write(message)

        self.logfile.write(message)

    def flush(self):
        if self.terminal:
            self.terminal.flush()

        if not self.logfile.closed:
            self.logfile.flush()

    def close(self):
        """
        Restore stdout, close current log file,
        and move all .log files from the current directory to log_dir.
        """
        if self._closed:
            return

        self._closed = True

        # Restore stdout first
        sys.stdout = self.terminal

        # Flush and close current log
        self.flush()

        if not self.logfile.closed:
            self.logfile.close()

        current_dir = os.getcwd()
        log_dir_abs = os.path.abspath(self.log_dir)

        for logfile in glob.glob(os.path.join(current_dir, "*.log")):
            src = os.path.abspath(logfile)

            # Skip files already inside log_dir
            if os.path.dirname(src) == log_dir_abs:
                continue

            dst = os.path.join(
                log_dir_abs,
                os.path.basename(src)
            )

            try:
                shutil.move(src, dst)
            except Exception as e:
                print(f"Could not move {src}: {e}")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
