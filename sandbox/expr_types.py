#
# Copyright 2014-2016 Vinay Vasista, Ravi Teja Mullapudi, Uday Bondhugula,
# and others from Multicore Computing Lab, Department of Computer Science
# and Automation, Indian Institute of Science
#

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#
# expr_types.py : Expression types and type precedence.
#

from __future__ import absolute_import, division, print_function

class Void(object):
    @staticmethod
    def c_type_name():
        return 'void'
    pass
class Float(object):
    @staticmethod
    def c_type_name():
        return 'float'
    pass
class Double(object):
    @staticmethod
    def c_type_name():
        return 'double'
    pass
class Int(object):
    @staticmethod
    def c_type_name():
        return 'int32'
    pass
class UInt(object):
    @staticmethod
    def c_type_name():
        return 'uint32'
    pass
class Short(object):
    @staticmethod
    def c_type_name():
        return 'int16'
    pass
class UShort(object):
    @staticmethod
    def c_type_name():
        return 'uint16'
    pass
class Char(object):
    @staticmethod
    def c_type_name():
        return 'int8'
    pass
class UChar(object):
    @staticmethod
    def c_type_name():
        return 'uint8'
    pass
class Long(object):
    @staticmethod
    def c_type_name():
        return 'int64'
    pass
class ULong(object):
    @staticmethod
    def c_type_name():
        return 'uint64'
    pass
class Rational(object):
    pass

def result_type(a, b):
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
