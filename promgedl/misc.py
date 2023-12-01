
import logging, sys

class colorcode():

    reset = "\033[0m"
    bold = "\033[1m"
    underline = "\033[4m"
    noUnderline = "\033[24m"
    reverseText = "\033[7m"
    noReverseText = "\033[27m"
    fgBlack = "\033[30m"
    bgBlack = "\033[40m"
    fgDarkRed = "\033[31m"
    bgDarkRed = "\033[41m"
    fgDarkGreen = "\033[32m"
    bgDarkGreen = "\033[42m"
    fgDarkYellow = "\033[33m"
    bgDarkYellow = "\033[43m"
    fgDarkBlue = "\033[34m"
    bgDarkBlue = "\033[44m"
    fgDarkMagenta = "\033[35m"
    bgDarkMagenta = "\033[45m"
    fgDarkCyan = "\033[36m"
    bgDarkCyan = "\033[46m"
    fgDarkWhite = "\033[37m"
    bgDarkWhite = "\033[47m"
    fgBrightBlack = "\033[90m"
    bgBrightBlack = "\033[100m"
    fgBrightRed = "\033[91m"
    bgBrightRed = "\033[101m"
    fgBrightGreen = "\033[92m"
    bgBrightGreen = "\033[102m"
    fgBrightYellow = "\033[93m"
    bgBrightYellow = "\033[103m"
    fgBrightBlue = "\033[94m"
    bgBrightBlue = "\033[104m"
    fgBrightMagenta = "\033[95m"
    bgBrightMagenta = "\033[105m"
    fgBrightCyan = "\033[96m"
    bgBrightCyan = "\033[106m"
    fgBrightWhite = "\033[97m"
    bgBrightWhite = "\033[107m"

class ColorFormatter(logging.Formatter):
    """Logging Formatter to add colors and count warning / errors"""
    formatINFO = " %(levelname)s - %(name)s - %(message)s"
    formatDEBUG = " %(levelname)s - (%(filename)s:%(lineno)d) - %(name)s :\t %(message)s "
    formatWARNING = " %(levelname)s - (%(filename)s:%(lineno)d) - %(name)s : %(message)s "
    formatERROR = " %(levelname)s - (%(filename)s:%(lineno)d) - %(name)s : %(message)s "
    formatCRITICAL = " %(levelname)s - (%(filename)s:%(lineno)d) %(asctime)s - %(name)s - %(message)s "


    FORMATS = {
        logging.DEBUG: colorcode.fgDarkCyan + colorcode.bold + formatDEBUG + colorcode.reset,
        logging.INFO: colorcode.fgBrightBlack + formatINFO + colorcode.reset,
        logging.WARNING: colorcode.fgBrightYellow + formatWARNING + colorcode.reset,
        logging.ERROR: colorcode.fgDarkRed + formatERROR + colorcode.reset,
        logging.CRITICAL: colorcode.bold + colorcode.fgBrightRed + formatCRITICAL + colorcode.reset
    }
    def __init__(self) -> None:
        super().__init__()

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)
    
def configureHandler():
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    ch.setFormatter(ColorFormatter())
    return ch

if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG, handlers=[configureHandler()])
    logger = logging.getLogger(__name__)
    logger.debug("Test Debug")
    logger.info("Test Info")
    logger.warning("Test Warning")
    logger.error("Test Error")
    logger.critical("Test Critical")