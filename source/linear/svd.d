module linear.svd;

import linear;
import std.math : sqrt, eq = approxEqual;
import std.numeric : dotProduct;
version (unittest) import std.stdio;

/** Singular value decomposition.

Params:
    A = the matrix to decompose.
Returns:
A struct with three members: U and tV of the same type as A and S of Diagonal!T.
U S tV = A holds.
*/
SVD!(T, rowMajor) svd(T, RowMajor rowMajor)(Matrix!(T, rowMajor) A)
{
    return SVD!(T, rowMajor)(A);
}

unittest
{
    [[1.0L, 3, 1], [3.0L, -2, -8], [-2.0L, 4, 8], [1.0L, 0, 1]].matrix.svd.original;
    [[1.0L, 3, 1], [3.0L, -2, -8], [-2.0L, 4, 8], [1.0L, 0, 1]].transpose.matrix.svd.original;
    [[1.0L, 3, 1], [3.0L, -2, -8], [-2.0L, 4, 8], [1.0L, 0, 1]].matrix!columnMajor.svd.original;
    [[1.0L, 3, 1], [3.0L, -2, -8], [-2.0L, 4, 8], [1.0L, 0, 1]].transpose.matrix!columnMajor.svd.original;
}

/// Cholesky-Banachiewicz decomposition of a symmetric matrix.
MatrixType choleskyBanachiewicz(MatrixType)(MatrixType A)
    if (is (MatrixType == Matrix!(T, rowMajor), T, RowMajor rowMajor))
in
{
    assert (A.rows == A.columns);
}
out (L)
{
    assert (A.approxEqual(L.timesAdjoint));
}
body
{
    immutable n = A.rows;
    auto L = MatrixType(n, n);
    foreach (i; 0..n)
    {
        foreach (j; 0..i)
            L[i, j] = (A[i, j] - L[i, 0..j] * L[j, 0..j].transpose) / L[j, j];
        L[i, i] = (A[i, i] - L[i, 0..i] * L[i, 0..i].transpose).sqrt;
    }
    return L;
}

/// Singular values of a matrix.
auto singularValues(MatrixType)(MatrixType A)
    if (is (MatrixType == Matrix!(T, rowMajor), T, RowMajor rowMajor))
{
    if (A.rows > A.columns)
        A = A.adjointTimes;
    else
        A = A.timesAdjoint;
    foreach (i; 0..10)
    {
        auto L = A.choleskyBanachiewicz;
        A = L.adjointTimes;
    }
    auto ret = A.diagonal;
    foreach (ref elem; ret.payload)
        elem = elem.sqrt;
    ret.payload.sort!((a, b) => a > b);
    return ret;
}

auto eigenvectors(MatrixType, U)(MatrixType A, U[] eigenvalues)
    if (is (MatrixType == Matrix!(T, rowMajor), T, RowMajor rowMajor) && is (typeof (MatrixType.init[0, 0] -= U.init)))
{
    return eigenvalues.map!(ev => A.eigenvectors(ev)).join.array;
}

/// Eigenvectors of a matrix.
auto eigenvectors(MatrixType, U)(MatrixType A, U eigenvalue)
    if (is (MatrixType == Matrix!(T, rowMajor), T, RowMajor rowMajor) && is (typeof (MatrixType.init[0, 0] -= U.init)))
in
{
    assert (A.rows == A.columns);
}
out (result)
{
    foreach (v; result)
        assert ((A * v.column).approxEqual(eigenvalue * v.column));
}
body
{
    immutable n = A.rows;
    static if (is (MatrixType == Matrix!(T, rowMajor), T, RowMajor rowMajor))
    {
    static if (rowMajor)
    auto a = A.copy.payload;
    else
    auto a = A.payload.transpose;
    }
    size_t[] dependent;
    foreach (i, row; a)
        row[i] -= eigenvalue;
    foreach (i; 0..n)
    {
        if (!a[i..$].findNonzeroRow)
            break;
        immutable ii = a[i].nonzeroIndex;
        dependent ~= ii;
        foreach (j; 0..n)
        {
            if (i == j)
                continue;
            immutable r = a[j][ii] / a[i][ii];
            a[j][] -= a[i][] * r;
        }
        immutable r = 1 / a[i][ii];
        a[i][] *= r;
    }
    auto free = n.iota.setDifference(dependent);
    MatrixType.Element[][] vectors;
    foreach (j; free)
    {
        auto v = new MatrixType.Element[n];
        v[] = 0;
        v[j] = 1;
        foreach (i, d; dependent)
            v[d] -= a[i][j];
        vectors ~= v;
    }
    return vectors;
}

/// Gram-Schmidt orthonormalization of vectors.
T[][] gramSchmidt(T)(T[][] a)
{
    immutable n = a.length;
    foreach (i; 0..n)
    {
        foreach (j; 0..i)
        {
            immutable r = a[i].dotProduct(a[j]) / a[j].dotProduct(a[j]);
            a[i][] -= a[j][] * r;
        }
        immutable r = 1 / a[i].dotProduct(a[i]).sqrt;
        a[i][] *= r;
    }
    return a;
}

private bool findNonzeroRow(T)(T[][] a)
{
    foreach (i; 0..a.length)
    {
        if (a[i].nonzero)
        {
            a[0].swap(a[i]);
            return true;
        }
    }
    return false;
}

private bool nonzero(T)(T[] a)
{
    foreach (elem; a)
        if (!elem.eq(0))
            return true;
    return false;
}

private size_t nonzeroIndex(T)(T[] a)
in
{
    assert (a.nonzero);
}
body
{
    foreach (i, elem; a)
        if (!elem.eq(0))
            return i;
    assert (false);
}

private:

struct SVD(T, RowMajor rowMajor)
{
    Matrix!(T, rowMajor) U, tV;
    Diagonal!T S;
    this (Matrix!(T, rowMajor) A)
    {
        S = A.singularValues;
        auto evalues = S.payload.map!(a => a * a).array;
        auto evectors = A.smallSquare.eigenvectors(evalues).gramSchmidt;
        if (A.rows <= A.columns)
        {
            U = evectors.transpose.matrix.to!rowMajor;
            tV = (A.deepTranspose * U / S).deepTranspose;
        }
        else
        {
            tV = evectors.matrix.to!rowMajor;
            U = (A * tV.deepTranspose) / S;
        }
    }
    Matrix!(T, rowMajor) original()
    {
        return U * S * tV;
    }
}

unittest
{
    import std.stdio;
    stderr.writeln("linear.svd: All green!");
}
