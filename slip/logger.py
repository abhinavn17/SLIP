import sys
import os
from datetime import datetime
import pytz

class PrintLogger:
    """
    Redirects all print statements to a log file and optionally to terminal.
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

        # Create log directory if it doesn't exist
        os.makedirs(log_dir, exist_ok=True)

        # Get current time in Asia/Kolkata timezone
        tz_IN = pytz.timezone('Asia/Kolkata')
        datetime_IN = datetime.now(tz_IN)
        dtobj = datetime_IN.strftime("%Y-%m-%d_%H%M%S")

        self.log_file = os.path.join(log_dir, f"{name}_{dtobj}.log")
        self.logfile = open(self.log_file, 'a', encoding='utf-8')

        # Terminal output
        self.terminal = sys.stdout if self.to_terminal else None

        sys.stdout = self

        self.write(f"Logging started. Log file: {self.log_file}\n")

    def write(self, message):
        if self.terminal:
            self.terminal.write(message)
        self.logfile.write(message)

    def flush(self):
        """Ensure compatibility with flush() calls."""
        if self.terminal:
            self.terminal.flush()
        self.logfile.flush()

    def close(self):
        """Close all log files."""
        self.logfile.close()