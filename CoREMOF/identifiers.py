"""CoRE release identifiers, distinct from chemical MOFid-v1/v2 strings.

Identifiers encode ``YYYY[elements][topology]dimension[variant]serial``.
The source database and access rights are metadata, not parts of a CoRE ID.
Published identifiers are persistent names: parsing does not reinterpret
their descriptive fields using newer calculations or infer group membership.
"""

from dataclasses import dataclass
import re


CORE_ID_SCHEMA_VERSION = "core-id/1.0"
CORE_ID_PATTERN = (
    r"(?P<year>[0-9]{4})\[(?P<elements>(?:[A-Z][a-z]?)+)\]"
    r"\[(?P<topology>[A-Za-z0-9][A-Za-z0-9_.-]*)\]"
    r"(?P<dimension>[0-3])\[(?P<variant>ASR|FSR|ION)\]"
    r"(?P<serial>[1-9][0-9]*)"
)
CORE_ID_RE = re.compile(CORE_ID_PATTERN)
_ELEMENTS = frozenset(
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe "
    "Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In "
    "Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf "
    "Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am "
    "Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og".split()
)
SOURCE_DATABASES = frozenset({"COD", "CSD", "SI"})


@dataclass(frozen=True)
class CoREID:
    """Parsed descriptive fields of an exact, persistent CoRE ID.

    ``year=0`` represents the literal unknown-year token ``0000``.
    ``dimension`` is bonded-framework dimensionality, not pore channels.
    ``elements`` includes metals or metalloids. ``topology='nan'`` is an
    identifier token, not an imputed scientific feature or a floating NaN.
    """

    year: int
    elements: str
    topology: str
    dimension: int
    variant: str
    serial: int

    def __str__(self):
        return format_core_id(
            self.year, self.elements, self.topology, self.dimension,
            self.variant, self.serial,
        )


def parse_core_id(value):
    """Parse an exact CoRE ID; reject paths, extensions and unsafe tokens."""
    match = CORE_ID_RE.fullmatch(value) if isinstance(value, str) else None
    if match is None:
        raise ValueError(
            "Use a CoRE ID in YYYY[elements][topology]dimension[variant]serial format"
        )
    fields = match.groupdict()
    symbols = re.findall(r"[A-Z][a-z]?", fields["elements"])
    if any(symbol not in _ELEMENTS for symbol in symbols) or len(set(symbols)) != len(symbols):
        raise ValueError("CoRE ID elements must be distinct valid element symbols")
    return CoREID(
        int(fields["year"]), fields["elements"], fields["topology"],
        int(fields["dimension"]), fields["variant"], int(fields["serial"]),
    )


def format_core_id(year, elements, topology, dimension, variant, serial):
    """Format validated fields without allocating serials or changing their order."""
    if type(year) is not int or not 0 <= year <= 9999:
        raise ValueError("year must be an integer in 0..9999 (0 means unknown)")
    if type(dimension) is not int or dimension not in (0, 1, 2, 3):
        raise ValueError("dimension must be an integer in 0..3")
    if type(serial) is not int or serial < 1:
        raise ValueError("serial must be a positive integer")
    if not all(isinstance(value, str) for value in (elements, topology, variant)):
        raise ValueError("elements, topology and variant must be exact strings")
    value = "{:04d}[{}][{}]{}[{}]{}".format(
        year, elements, topology, dimension, variant, serial
    )
    parse_core_id(value)
    return value


def validate_source_database(value):
    """Require an explicit source field; CoRE IDs do not encode a source."""
    if value not in SOURCE_DATABASES:
        raise ValueError("source_database must be explicit metadata: COD, CSD or SI")
    return value
