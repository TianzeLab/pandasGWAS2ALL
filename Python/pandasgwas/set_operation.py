from pandasgwas.Ancestry import Ancestry
from pandasgwas.Association import Association
from pandasgwas.Gene import Gene
from pandasgwas.GenomicContext import GenomicContext
from pandasgwas.GenomicContext import GenomicContext
from pandasgwas.Locus import Locus
from pandasgwas.Publication import Publication
from pandasgwas.Snp import Snp
from pandasgwas.Study import Study
from pandasgwas.Trait import Trait


def bind(
        a: Ancestry | Association | Gene | GenomicContext | GenomicContext | Locus | Publication | Snp | Study | Trait,
        b: Ancestry | Association | Gene | GenomicContext | GenomicContext | Locus | Publication | Snp | Study | Trait) ->Ancestry | Association | Gene | GenomicContext | GenomicContext | Locus | Publication | Snp | Study | Trait :
    """
    Binds together GWAS objects of the same object. Note that bind() preserves duplicates whereas union() does not.

    Args:
        a: An object of the pandasgwas custom class.
        b: An object of the same type as a.

    Returns:
        An object of the same type as a.
    ```python
    from pandasgwas.get_snps import get_snps
    from pandasgwas.set_operation import bind


    a = get_snps(rs_id='rs3093017')
    b = get_snps(rs_id='rs7671745')
    c = bind(a,b)
    ```

    """
    return a + b


def intersect(
        a: Ancestry | Association | Gene | GenomicContext | GenomicContext | Locus | Publication | Snp | Study | Trait,
        b: Ancestry | Association | Gene | GenomicContext | GenomicContext | Locus | Publication | Snp | Study | Trait) ->Ancestry | Association | Gene | GenomicContext | GenomicContext | Locus | Publication | Snp | Study | Trait :
    """
    Returns the data common to both A and B, with no repetitions

    Args:
        a: An object of the pandasGWAS custom class.
        b:An object of the same type as a.

    Returns:
        An object of the same type as a.
    ```python
    from pandasgwas.get_snps import get_snps
    from pandasgwas.set_operation import bind


    a = get_snps(rs_id='rs3093017')
    b = get_snps(rs_id='rs7671745')
    c = intersect(a,b)
    ```

    """
    return a & b


def set_diff(
        a: Ancestry | Association | Gene | GenomicContext | GenomicContext | Locus | Publication | Snp | Study | Trait,
        b: Ancestry | Association | Gene | GenomicContext | GenomicContext | Locus | Publication | Snp | Study | Trait) -> Ancestry | Association | Gene | GenomicContext | GenomicContext | Locus | Publication | Snp | Study | Trait:
    """
    returns the data in A that is not in B, with no repetitions

    Args:
        a: An object of the pandasGWAS custom class.
        b:An object of the same type as a.

    Returns:
        An object of the same type as a.
    ```python
    from pandasgwas.get_snps import get_snps
    from pandasgwas.set_operation import bind


    a = get_snps(rs_id='rs3093017')
    b = get_snps(rs_id='rs7671745')
    c = set_diff(a,b)
    ```

    """
    return a - b


def set_xor(
        a: Ancestry | Association | Gene | GenomicContext | GenomicContext | Locus | Publication | Snp | Study | Trait,
        b: Ancestry | Association | Gene | GenomicContext | GenomicContext | Locus | Publication | Snp | Study | Trait) -> Ancestry | Association | Gene | GenomicContext | GenomicContext | Locus | Publication | Snp | Study | Trait:
    """
    returns the data of A and B that are not in their intersection (the symmetric difference), with no repetitions

    Args:
        a: An object of the pandasGWAS custom class.
        b:An object of the same type as a.

    Returns:
        An object of the same type as a.
    ```python
    from pandasgwas.get_snps import get_snps
    from pandasgwas.set_operation import bind


    a = get_snps(rs_id='rs3093017')
    b = get_snps(rs_id='rs7671745')
    c = set_xor(a,b)
    ```

    """
    return a ^ b


def union(
        a: Ancestry | Association | Gene | GenomicContext | GenomicContext | Locus | Publication | Snp | Study | Trait,
        b: Ancestry | Association | Gene | GenomicContext | GenomicContext | Locus | Publication | Snp | Study | Trait) -> Ancestry | Association | Gene | GenomicContext | GenomicContext | Locus | Publication | Snp | Study | Trait:
    """
    returns the combined data from A and B with no repetitions

    Args:
        a: An object of the pandasGWAS custom class.
        b:An object of the same type as a.

    Returns:
        An object of the same type as a.
    ```python
    from pandasgwas.get_snps import get_snps
    from pandasgwas.set_operation import bind


    a = get_snps(rs_id='rs3093017')
    b = get_snps(rs_id='rs7671745')
    c = union(a,b)
    ```

    """
    return a | b


def set_equal(
        a: Ancestry | Association | Gene | GenomicContext | GenomicContext | Locus | Publication | Snp | Study | Trait,
        b: Ancestry | Association | Gene | GenomicContext | GenomicContext | Locus | Publication | Snp | Study | Trait) -> Ancestry | Association | Gene | GenomicContext | GenomicContext | Locus | Publication | Snp | Study | Trait:
    """
    Check if the raw data of a and b are equal

    Args:
        a: An object of the pandasGWAS custom class.
        b:An object of the same type as a.

    Returns:
        True, if a is equal to b. False, if a is not equal to b
    ```python
    from pandasgwas.get_snps import get_snps
    from pandasgwas.set_operation import bind


    a = get_snps(rs_id='rs3093017')
    b = get_snps(rs_id='rs7671745')
    c = set_equal(a,b)
    ```

    """

    return a == b
