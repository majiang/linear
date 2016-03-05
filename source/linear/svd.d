module linear.svd;

import linear;
import std.math : sqrt, abs, isNaN, hypot, sgn, eq = approxEqual;
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
    auto arr = [[1.0L, 3, 1], [3.0L, -2, -8], [-2.0L, 4, 8], [1.0L, 0, 1]];
    arr.matrix.svd.original;
    arr.transpose.matrix.svd.original;
    arr.matrix!columnMajor.svd.original;
    arr.transpose.matrix!columnMajor.svd.original;
    arr = [[2.0L, 0, 0], [0.0L, 2, 0], [0.0L, 0, 1], [0.0L, 0, 0]];
    arr.matrix.svd.original;
    arr.transpose.matrix.svd.original;
    arr.matrix!columnMajor.svd.original;
    arr.transpose.matrix!columnMajor.svd.original;
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

/// ditto
auto singularValuesJacobi(MatrixType)(MatrixType A, real eps = 1e-20)
    if (is (MatrixType == Matrix!(T, rowMajor), T, RowMajor rowMajor))
{
    auto ret = A.smallSquare.jacobi!MatrixType(eps);
    foreach (ref elem; ret.payload)
        elem = elem.sqrt;
    ret.payload.sort!((a, b) => a > b);
    return ret;
}


// L<sup><var>p</var></sup>-distance of two diagonal matrices.
private real distanceL(real p, T)(Diagonal!T lhs, Diagonal!T rhs)
{
    auto ret = real(0);
    assert (lhs.size == rhs.size);
    foreach (x; lhs.payload.zip(rhs.payload))
        ret += (x[0] - x[1]).abs ^^ p;
    return ret ^^ (1 / p);
}

/// Jacobi eigenvalue algorithm for symmetric matrix.
auto jacobi(MatrixType)(MatrixType A, MatrixType.Element eps = 1e-20)
{
    size_t p, q;
    while (eps < A.largest!MatrixType(p, q))
        Rotation!(MatrixType.Element)(A.rows, p, q, A[p, p], A[p, q], A[q, q]).apply(A);
    return A.diagonal;
}

private M.Element largest(M)(M A, out size_t p, out size_t q)
{
    M.Element elem = 0;
    foreach (i; 0..A.rows)
        foreach (j; 0..A.columns)
            if (i != j && elem < A[i, j].abs)
            {
                p = i;
                q = j;
                elem = A[i, j].abs;
            }
    return elem;
}

private struct Rotation(T)
{
    immutable size_t size, p, q;
    immutable T app, apq, aqq, c, s;
    invariant
    {
        assert (!c.isNaN);
        assert (!s.isNaN);
    }
    this (in size_t size, in size_t p, in size_t q, in T app, in T apq, in T aqq)
    in
    {
        assert (p != q);
        assert (p < size);
        assert (q < size);
        assert (!app.isNaN);
        assert (!apq.isNaN);
        assert (!aqq.isNaN);
        assert (apq != 0);
    }
    body
    {
        this.size = size;
        this.p = p;
        this.q = q;
        this.app = app;
        this.apq = apq;
        this.aqq = aqq;
        immutable
            alpha = (app - aqq) / 2,
            beta = -apq,
            gamma = alpha.abs / alpha.hypot(beta);
        assert (!alpha.isNaN);
        assert (! beta.isNaN);
        assert (!gamma.isNaN);
        this.c = ((1 + gamma) / 2).sqrt;
        this.s = ((1 - gamma) / 2).sqrt * (alpha * beta).sgn;
    }
    void apply(MatrixType)(MatrixType A)
    {
        foreach (i; 0..size)
        {
            immutable tmp = c * A[p, i] - s * A[q, i];
            A[q, i] = s * A[p, i] + c * A[q, i];
            A[p, i] = tmp;
        }
        foreach (i; 0..size)
        {
            A[i, p] = A[p, i];
            A[i, q] = A[q, i];
        }
        immutable crossTerm = 2 * c * s * apq;
        A[p, p] = c * c * app + s * s * aqq - crossTerm;
        A[q, q] = s * s * app + c * c * aqq + crossTerm;
        A[p, q] = 0;
        A[q, p] = 0;
    }
}

/** Eigenvectors of a matrix.

Params:
    A = the matrix.
    eigenvalues = the eigenvalues of a matrix, multiplicity must be spcified by duplicating the values.
*/
MatrixType.Element[][] eigenvectors(MatrixType, U)(MatrixType A, U[] eigenvalues)
    if (is (MatrixType == Matrix!(T, rowMajor), T, RowMajor rowMajor) && is (typeof (MatrixType.init[0, 0] -= U.init)))
out (result)
{
    debug (evec)
        tracef("eigenvectors found:\n%(  %(%e %)\n%)\n.", result);
}
body
{
    MatrixType.Element[][] ret;
    size_t skip;
    foreach (ev; eigenvalues)
    {
        if (skip)
        {
            trace("skip");
            skip -= 1;
            continue;
        }
        auto tmp = A.eigenvectors(ev);
        skip = tmp.length - 1;
        ret ~= tmp;
    }
    return ret;
}

/// Eigenvectors of a matrix for a given eigenvalue.
MatrixType.Element[][] eigenvectors(MatrixType, U)(MatrixType A, U eigenvalue)
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
    debug (evec)
        tracef("eigenvectors found for eigenvalue %e:\n%(  %(%e %)\n%)\n.", eigenvalue, result);
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
in
{
    debug (gson)
        tracef("gson input:\n%(  %(%e %)\n%)\n.", a);
}
out (result)
{
    foreach (i, u; result)
        foreach (j, v; result[0..i+1])
            assert (u.dotProduct(v).eq((i == j)),
                "gson: %d-th and %d-th vector of output has product %f\n  %(%e %)\nand\n  %(%e %)\n.".format(
                    i, j, u.dotProduct(v), u, v));
}
body
{
    immutable n = a.length;
    foreach (i; 0..n)
    {
        foreach (j; 0..i)
        {
            debug (gson) tracef("gson: subtract %d-th vector [%(%e %)] from %d-th vector [%(%e %)]", j, a[j], i, a[i]);
            immutable r = a[i].dotProduct(a[j]) / a[j].dotProduct(a[j]);
            a[i][] -= a[j][] * r;
        debug (gson)
            tracef("gson now:\n%(  %(%e %)\n%)\n.", a);
        }
        auto rinv = a[i].dotProduct(a[i]);
        criticalf(rinv == 0, "gson: %d-th vector is now zero", i);
        rinv = rinv.sqrt;
        a[i][] /= rinv;
        debug (gson)
            tracef("gson now:\n%(  %(%e %)\n%)\n.", a);
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
    out
    {
        assert (original.approxEqual(A), "SVD equality does not hold:\n%(  %(%e %)\n%)\nand\n%(  %(%e %)\n%)\n.%s".format(original, A, this));
    }
    body
    {
        S = A.singularValuesJacobi;
        auto evalues = S.payload.map!(a => a * a).array;
        debug (svd)
            tracef("evalues: %(%e %)", evalues);
        auto evectors = A.smallSquare.eigenvectors(evalues).gramSchmidt;
        debug (svd)
            tracef("evectors:\n%(  %(%e %)\n%)\n.", evectors);
        if (A.rows < A.columns)
        {
            debug (svd)
                trace("rows <= columns");
            U = evectors.transpose.matrix.to!rowMajor;
            tV = (A.deepTranspose * U / S).deepTranspose;
        }
        else
        {
            debug (svd)
                trace("columns <= rows");
            tV = evectors.matrix.to!rowMajor;
            U = (A * tV.deepTranspose) / S;
        }
    }
    Matrix!(T, rowMajor) original()
    {
        return U * S * tV;
    }
    string toString()
    {
        return "SVD
(
U = %(%(%e %)\n    %)\n
S = %s\n
tV= %(%(%e %)\n    %)\n
)".format(U, S, tV);
    }
}

unittest
{
    import std.stdio;
    stderr.writeln("linear.svd: All green!");
}
