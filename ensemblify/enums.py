import enum


class StringEnum(enum.Enum):
    def __str__(self):
        return self.value

    def __add__(self, other):
        return str(self) + other

    def __radd__(self, other):
        return other + str(self)

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return str(self) == str(other)

    def __hash__(self):
        return hash(str(self))


# # # # # # # Data Types # # # # # # #


# class DataTypeEnum(StringEnum):
#     pass
#     DNA = "DNA"
#     AA = "AA"



