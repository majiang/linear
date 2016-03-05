module linear.svd;

import linear;
import std.math : sqrt, abs, isNaN, eq = approxEqual;
import std.numeric : dotProduct;
debug (sval) import std.stdio;
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
        {
            if (L[j, j] == 0)
            {
                debug (cbd)
                    tracef("cbd: division by L[%d, %d] = 0", j, j);
                immutable tmp = A[i, j] - L[i, 0..j] * L[j, 0..j].transpose;
                if (tmp != 0)
                    warning("cbd: treat L[%d, %d] = %e / 0 as 0", i, j, tmp);
                continue;
            }
            L[i, j] = (A[i, j] - L[i, 0..j] * L[j, 0..j].transpose) / L[j, j];
        }
        immutable tmp = A[i, i] - L[i, 0..i] * L[i, 0..i].transpose;
        if (tmp < 0)
        {
            warning("cbd: treat L[%d, %d] = sqrt(%e) as 0", i, i, tmp);
            continue;
        }
        L[i, i] = tmp.sqrt;
    }
    return L;
}

/// Singular values of a matrix.
auto singularValuesQR(MatrixType)(MatrixType A, real eps = 1e-20)
    if (is (MatrixType == Matrix!(T, rowMajor), T, RowMajor rowMajor))
{
    A = A.smallSquare;
    auto diff = MatrixType.Element.infinity;
    while (eps < diff)
    {
        auto L = A.choleskyBanachiewicz;
        auto old = A.diagonal;
        A = L.adjointTimes;
        diff = A.diagonal.distanceL!1(old);
        debug (svqr)
            tracef("error = %e", diff);
    }
    auto ret = A.diagonal;
    foreach (ref elem; ret.payload)
        elem = elem.sqrt;
    ret.payload.sort!((a, b) => a > b);
    return ret;
}

private real distanceL(real p, T)(Diagonal!T lhs, Diagonal!T rhs)
{
    auto ret = real(0);
    assert (lhs.size == rhs.size);
    foreach (x; lhs.payload.zip(rhs.payload))
        ret += (x[0] - x[1]).abs ^^ p;
    return ret ^^ (1 / p);
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
    assert (result.length,
        "No eigenvector found for eigenvalue %e of matrix:\n%(  %(%e %)\n%)\n.".format(
            eigenvalue, A));
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
            debug (gson) tracef("gson: subtract %d-th vector [%(%e %)] from %d-th vector [%(%e %)]", j, a[j], i, a[j]);
            immutable r = a[i].dotProduct(a[j]) / a[j].dotProduct(a[j]);
            a[i][] -= a[j][] * r;
        }
        immutable rinv = a[i].dotProduct(a[i]);
        criticalf(rinv == 0, "gson: %d-th vector is now zero", i);
        a[i][] /= rinv;
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
        S = A.singularValuesQR;
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
