"""Fragment canonicalisation: attachment-point rooting and multi-dummy forms."""
from isogrow.bioisostere.mmpdb_backend import canonical_fragment, fragment_forms


def test_benzofuran_is_rooted_at_the_dummy():
    # mmpdb stores the dummy-rooted form; a plain RDKit canonicalisation would
    # reorder this to c1ccc2oc([*:1])cc2c1 and never match.
    assert canonical_fragment("[*:1]c1cc2ccccc2o1") == "[*:1]c1cc2ccccc2o1"
    assert canonical_fragment("c1ccc2oc([*:1])cc2c1") == "[*:1]c1cc2ccccc2o1"


def test_phenyl_and_heteroaryl():
    assert canonical_fragment("[*:1]c1ccccc1") == "[*:1]c1ccccc1"
    assert canonical_fragment("c1ccc([*:1])cc1") == "[*:1]c1ccccc1"


def test_invalid_smiles_returns_none():
    assert canonical_fragment("not-a-smiles") is None
    assert canonical_fragment("") is None


def test_fragment_forms_cover_all_roots():
    single = fragment_forms("[*:1]c1ccc(OC)cc1")
    assert "[*:1]c1ccc(OC)cc1" in single

    linker = "[*:1]C(=O)Oc1ccc([*:2])cc1"
    forms = fragment_forms(linker)
    # raw + plain canonical + one rooted form per attachment point
    assert linker in forms
    assert len(forms) >= 3


def test_fragment_forms_empty_input():
    assert fragment_forms("") == set()
