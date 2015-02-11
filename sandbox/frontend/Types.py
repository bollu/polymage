# Making things compatible for python 3
# Yet to figure out how to make range like xrange
from __future__ import absolute_import, division, print_function

class Void(object):
    @staticmethod
    def cTypeName():
        return 'void'
    pass
class Float(object):
    @staticmethod
    def cTypeName():
        return 'float'
    pass
class Double(object):
    @staticmethod
    def cTypeName():
        return 'double'
    pass
class Int(object):
    @staticmethod
    def cTypeName():
        return 'int32'
    pass
class UInt(object):
    @staticmethod
    def cTypeName():
        return 'uint32'
    pass
class Short(object):
    @staticmethod
    def cTypeName():
        return 'int16'
    pass
class UShort(object):
    @staticmethod
    def cTypeName():
        return 'uint16'
    pass
class Char(object):
    @staticmethod
    def cTypeName():
        return 'int8'
    pass
class UChar(object):
    @staticmethod
    def cTypeName():
        return 'uint8'
    pass
class Long(object):
    @staticmethod
    def cTypeName():
        return 'int64'
    pass
class ULong(object):
    @staticmethod
    def cTypeName():
        return 'uint64'
    pass
class Rational(object):
    pass

def implictTypeConversion(a, b):
    if a is Double or b is Double:
        return Double
    elif a is Float or b is Float:
        return Float
    elif a is ULong or b is ULong:
        return ULong
    elif a is Long or b is Long:
        return Long
    elif a is UInt or b is UInt:
        return UInt
    elif a is Int or b is Int:
        return Int
    elif a is UShort or b is UShort:
        return UShort
    elif a is Short or b is Short:
        return Short
    elif a is UChar or b is UChar:
        return UChar
    elif a is Char and b is Char:
        return Char
    raise TypeError((a, b))
