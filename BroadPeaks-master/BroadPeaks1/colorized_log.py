#!/usr/bin/env/ python
import logging

class ColoredFormatter(logging.Formatter):
    """Special custom formatter for colorizing log messages!"""

    BLACK = '\033[0;30m'
    RED = '\033[0;31m'
    GREEN = '\033[0;32m'
    BROWN = '\033[0;33m'
    BLUE = '\033[0;34m'
    PURPLE = '\033[0;35m'
    CYAN = '\033[0;36m'
    GREY = '\033[0;37m'

    DARK_GREY = '\033[1;30m'
    LIGHT_RED = '\033[1;31m'
    LIGHT_GREEN = '\033[1;32m'
    YELLOW = '\033[1;33m'
    LIGHT_BLUE = '\033[1;34m'
    LIGHT_PURPLE = '\033[1;35m'
    LIGHT_CYAN = '\033[1;36m'
    WHITE = '\033[1;37m'

    RESET = "\033[0m"

    def __init__(self, *args, **kwargs):
        self._colors = {logging.DEBUG: self.DARK_GREY,
                        logging.INFO: self.BROWN,
                        logging.WARNING: self.GREEN,
                        logging.ERROR: self.RED,
                        logging.CRITICAL: self.LIGHT_RED}
        super(ColoredFormatter, self).__init__(*args, **kwargs)

    def format(self, record):
        """Applies the color formats"""
        record.msg = self._colors[record.levelno] + record.msg + self.RESET
        return logging.Formatter.format(self, record)

    def setLevelColor(self, logging_level, escaped_ansi_code):
        self._colors[logging_level] = escaped_ansi_code

"""
console = logging.StreamHandler()
console.setLevel(logging.DEBUG)
formatter = ColoredFormatter('%(message)s')
console.setFormatter(formatter)
logging.getLogger('').addHandler(console)

# messages below WARN do not supported why???

logging.warn('df')
logging.error('err')
logging.info('drtwesaf')
logging.debug('sf')
logging.critical('sdf')
"""
