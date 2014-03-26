from __future__ import print_function
import logging


class EnrichError(Exception):
    """
    Exception class for errors in the Enrich2 project.
    """
    def __init__(self, value, name):
        self.value = value
        self.name = name
        logging.error(str(self))

    def __str__(self):
        return "{message!r} [{name}]".format(message=self.value, name=self.name)



class EnrichMessage(object):
    """
    Class for warning messages in the Enrich2 project.

    .. warning:: Not currently used.
    """
    def __init__(self, message):
        self.message = message
    def __str__(self):
        return repr(self.message)
