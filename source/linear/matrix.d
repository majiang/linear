module linear.matrix;

import linear;
import std.typecons;
import std.traits : Unqual;
alias RowMajor = Flag!"RowMajor";
enum rowMajor = RowMajor.yes, columnMajor = RowMajor.no;

/** Matrix type. The layout can be chosen from row-major and column-major.
*/
struct Matrix(T, RowMajor rowMajor=rowMajor)
    if (is (T == Unqual!T))
{
    /// Initialize Matrix with given payload.
    this (T[][] payload)
    {
        this.payload = payload;
    }
    /// Initialize zero Matrix with specified rows and columns;
    static if (rowMajor)
    this (size_t rows, size_t columns)
    {
        auto _payload = new T[rows * columns];
        static if (T.init != 0)
            _payload[] = 0;
        foreach (i; 0..rows)
            payload ~= _payload[columns * i .. columns * (i+1)];
    }
    /// ditto
    static if (!rowMajor)
    this (size_t rows, size_t columns)
    {
        auto _payload = new T[rows * columns];
        static if (T.init != 0)
            _payload[] = 0;
        foreach (i; 0..columns)
            payload ~= _payload[rows * i .. rows * (i+1)];
    }
    /// Addition and subtraction.
    auto opOpAssign(string op)(in typeof (this) rhs)
        if (op == "+" || op == "-")
    {
        foreach (i, major; rhs.payload)
    mixin ("payload[i][] "~op~"= major[];");
        return this;
    }
    /// ditto
    auto opBinary(string op)(in typeof (this) rhs) const
        if (op == "+" || op == "-")
    {
        T[][] payload;
        foreach (major; this.payload)
            payload ~= major.dup;
mixin ("return (Unqual!(typeof (this))(payload)) "~op~"= rhs;");
    }
    /// Scaler multiplication.
    auto opOpAssign(string op)(in T rhs)
        if (op == "*" || op == "/")
    {
        foreach (major; payload)
    mixin ("major[] "~op~"= rhs;");
        return this;
    }
    /// ditto
    auto opBinary(string op)(in T rhs) const
        if (op == "*" || op == "/")
    {
        return copy.opOpAssign!op(rhs);
    }
    /// ditto
    auto opBinaryRight(string op)(in T lhs) const
        if (op == "*")
    {
        return opBinary!op(lhs);
    }
    /// Multiplication with vector.
    auto opBinaryRight(string op)(in RowVector!T lhs) const
        if (op == "*" && !rowMajor)
    {
        RowVector!T ret;
        foreach (column; this.payload)
            ret.payload ~= lhs * column;
        return ret;
    }
    /// ditto
    auto opBinaryRight(string op)(in RowVector!T lhs) const
        if (op == "*" && rowMajor)
    in
    {
        assert (payload.length == lhs.payload.length);
    }
    body
    {
        auto ret = RowVector!T(this.minorLength());
        foreach (i, row; this.payload)
            ret.payload[] += lhs.payload[i] * row[];
        return ret;
    }
    /// ditto
    auto opBinary(string op)(in ColumnVector!T rhs) const
        if (op == "*" && rowMajor)
    {
        ColumnVector!T ret;
        foreach (row; this.payload)
            ret.payload ~= row * rhs;
        return ret;
    }
    /// ditto
    auto opBinary(string op)(in ColumnVector!T rhs) const
        if (op == "*" && !rowMajor)
    in
    {
        assert (payload.length == rhs.payload.length);
    }
    body
    {
        auto ret = ColumnVector!T(this.minorLength());
        static if (T.init != 0)
            ret.payload[] = 0;
        foreach (i, column; this.payload)
            ret.payload[] += column[] * rhs.payload[i];
        return ret;
    }
    /// Matrix multiplication.
    auto opBinary(string op)(in typeof (this) rhs) const
        if (op == "*" && rowMajor)
    {
        Unqual!(typeof (this)) ret;
        foreach (row; payload)
            ret.payload ~= (const RowVector!T(row) * rhs).payload;
        return ret;
    }
    /// ditto
    auto opBinary(string op)(in typeof (this) rhs) const
        if (op == "*" && !rowMajor)
    {
        Unqual!(typeof (this)) ret;
        foreach (column; rhs.payload)
            ret.payload ~= (this * const ColumnVector!T(column)).payload;
        return ret;
    }
    /// The number of rows.
    auto rows() const
    {
        static if (rowMajor)
            return majorLength;
        else
            return minorLength;
    }
    /// The number of columns.
    auto columns() const
    {
        static if (rowMajor)
            return minorLength;
        else
            return majorLength;
    }
    /// dollar operator.
    auto opDollar(size_t dimension)() const
        if (dimension < 2)
    {
        static if (dimension == 0)
            return rows;
        else
            return columns;
    }
    Tuple!(size_t, size_t) opSlice(size_t dimension)(in size_t min, in size_t max) const
        if (dimension < 2)
    in
    {
        assert (min <= max);
        static if (dimension == 0)
            assert (max <= rows);
        else
            assert (max <= columns);
    }
    body
    {
        return min.tuple!(size_t, size_t)(max);
    }
    RowVector!T opIndex(size_t i, Tuple!(size_t, size_t) jrange) const
    {
        static if (rowMajor)
            return row(i)[jrange[0]..jrange[1]];
        else
            return getRow(i)[jrange[0]..jrange[1]];
    }
    ColumnVector!T opIndex(Tuple!(size_t, size_t) irange, size_t j) const
    {
        static if (!rowMajor)
            return column(j)[irange[0]..irange[1]];
        else
            return getColumn(j)[irange[0]..irange[1]];
    }
    Matrix!(T, rowMajor) opIndex(Tuple!(size_t, size_t) irange, Tuple!(size_t, size_t) jrange) const
    {
        static if (rowMajor)
            return payload.minor(irange[0], irange[1], jrange[0], jrange[1]).matrix;
        else
            return payload.minor(jrange[0], jrange[1], irange[0], irange[1]).matrix!columnMajor;
    }
    /// Get a row/column.
    static if (rowMajor)
    auto row(in size_t index) inout
    {
        return payload[index].row;
    }
    /// ditto
    static if (!rowMajor)
    auto column(in size_t index) inout
    {
        return payload[index].column;
    }
    /// Iterate over rows/columns.
    static if (rowMajor)
    auto byRow()
    {
        struct R
        {
            T[][] payload;
            auto front()
            {
                return payload.front.row;
            }
            void popFront()
            {
                payload.popFront;
            }
            auto empty()
            {
                return payload.empty;
            }
        }
        return R(payload);
    }
    /// ditto
    static if (!rowMajor)
    auto byColumn()
    {
        struct R
        {
            T[][] payload;
            auto front()
            {
                return payload.front.column;
            }
            void popFront()
            {
                payload.popFront;
            }
            auto empty()
            {
                return payload.empty;
            }
        }
        return R(payload);
    }
    /// Get a row/column.
    static if (rowMajor)
    auto getColumn(in size_t index) const
    {
        T[] ret;
        foreach (row; payload)
            ret ~= row[index];
        return ret.column;
    }
    /// ditto
    static if (!rowMajor)
    auto getRow(in size_t index) const
    {
        T[] ret;
        foreach (column; payload)
            ret ~= column[index];
        return ret.row;
    }
    /// O(1) transposition.
    auto transpose()
    {
        return Matrix!(T, cast(RowMajor)(!rowMajor))(payload);
    }
    /// O(rc) transposition.
    auto deepTranspose() const
    {
        return Matrix!(T, rowMajor)(payload.transpose);
    }
    /// deep copy.
    auto copy() const
    {
        T[][] ret;
        foreach (major; payload)
            ret ~= major.dup;
        return Unqual!(typeof (this))(ret);
    }
    /// For element-wise special operation.
    auto _rawPayload()
    {
        return payload;
    }
    /// ditto
    auto ref inout(T) opIndex(size_t i, size_t j) inout
    {
        static if (rowMajor)
            return payload[i][j];
        else
            return payload[j][i];
    }
    import std.format;
    void toString(scope void delegate(const(char)[]) sink, FormatSpec!char fmt) const
    {
        sink.formatValue(rowMajor ? payload : payload.transpose, fmt);
    }
package:
    T[][] payload;
private:
    size_t majorLength() const
    {
        return payload.length;
    }
    size_t minorLength() const
    {
        return majorLength ? payload.front.length : 0;
    }
}

/// Construct matrix from payload.
auto matrix(RowMajor rowMajor=rowMajor, T)(T[][] payload)
{
    return Matrix!(T, rowMajor)(payload);
}

unittest
{
    auto a = Matrix!double(1, 1);
    assert (a.payload == [[0.0]]);
}
unittest
{
    auto a = Matrix!(int, rowMajor)(2, 3);
    a.payload[1][] = 3;
    auto b = [[0, 1, 2], [0, 1, 2]].matrix;
    assert (a.rows == b.rows);
    assert (a.columns == b.columns);
    assert ((a + b).payload == [[0, 1, 2], [3, 4, 5]]);
    assert ((a - b).payload == [[0,-1,-2], [3, 2, 1]]);
}
unittest
{
    auto a = Matrix!(int, columnMajor)(3, 2);
    a.payload[1][] = 3;
    auto b = [[0, 1, 2], [0, 1, 2]].matrix!columnMajor;
    assert (a.rows == b.rows);
    assert (a.columns == b.columns);
    assert ((a + b).payload == [[0, 1, 2], [3, 4, 5]]);
    assert ((a - b).payload == [[0,-1,-2], [3, 2, 1]]);
}
debug (memoryUsage) unittest
{
    Matrix!int[100] x;
    foreach (i; 0..100)
        x[i] = Matrix!int(100000, 10);
    import std.stdio;
    stderr.writeln("...");
    readln;
    // 724388KB -> 490412KB
}
unittest
{
    auto a = [[0, 1, 2], [3, 4, 5]].matrix;
    auto b = [0, 1, 2].column;
    assert ((a * b).payload == [5, 14]);
}
unittest
{
    auto a = [0, 1, 2].row;
    auto b = [[0, 1, 2], [3, 4, 5]].matrix!columnMajor;
    assert ((a * b).payload == [5, 14]);
}
unittest
{
    auto a = [0, 1, 2].row;
    auto b = [[0, 1], [2, 3], [4, 5]].matrix;
    assert ((a * b).payload == [10, 13]);
}
unittest
{
    auto a = [[0, 1], [2, 3], [4, 5]].matrix!columnMajor;
    auto b = [0, 1, 2].column;
    assert ((a * b).payload == [10, 13]);
}
unittest
{
    auto a = [[0, 1, 2], [3, 4, 5]].matrix;
    auto b = [[6, 7], [8, 9], [10, 11]].matrix;
    assert ((a * b).payload == [[28, 31], [100, 112]]);
}
unittest
{
    auto a = [[0, 1, 2], [3, 4, 5]].matrix!columnMajor;
    auto b = [[6, 7], [8, 9], [10, 11]].matrix!columnMajor;
    assert ((a * b).payload == [[21, 34, 47], [27, 44, 61], [33, 54, 75]]);
}
unittest
{
    import std.exception, core.exception;
    auto a = [[0, 1], [2, 3]].matrix;
    auto b = [[0, 1, 2], [3, 4, 5]].matrix;
    auto c = [[0, 1, 2], [3, 4, 5], [6, 7, 8]].matrix;
    assert ((a * b * c).payload == [[42, 54, 66], [156, 198, 240]]);
    assertThrown!AssertError(b * a);
    assertThrown!AssertError(c * b);
    auto p = c.transpose;
    auto q = b.transpose;
    auto r = a.transpose;
    assert ((p * q * r).payload == [[42, 54, 66], [156, 198, 240]]);
    assertThrown!AssertError(q * p);
    assertThrown!AssertError(r * q);
    a = a.deepTranspose;
    b = b.deepTranspose;
    c = c.deepTranspose;
    assert ((c * b * a).payload == [[42, 156], [54, 198], [66, 240]]);
    assertThrown!AssertError(a * b);
    assertThrown!AssertError(b * c);
}
unittest
{
    auto a = [[0, 0], [0, 0]].matrix;
    auto b = [[0, 0], [0, 0]].matrix!columnMajor;
    a[0, 1] = 1;
    b[0, 1] = 1;
    assert (a.payload[0][1] == 1);
    assert (b.payload[1][0] == 1);
    const p = a;
    const q = b;
    assert (p[0, 1] == 1);
    assert (q[0, 1] == 1);
}
unittest
{
    auto a = [[0, 1, 2], [3, 4, 5], [6, 7, 8]].matrix;
    assert (a[0, 0] == 0);
    assert (a[0, 1..$] == [1, 2].row);
    assert (a[1..$, 1..$] == [[4, 5], [7, 8]].matrix);
    assert (a[1..$, 0] == [3, 6].column);
}

/// transpose a payload for matrix.
T[][] transpose(T)(in T[][] payload)
{
    if (payload.empty || payload.front.empty)
        return [];
    auto ret = new T[][](payload[0].length, payload.length);
    foreach (i, major; payload)
        foreach (j, elem; major)
            ret[j][i] = elem;
    return ret;
}

unittest
{
    assert([[0, 1], [2, 3]].transpose == [[0, 2], [1, 3]]);
}

unittest
{
    auto a = Matrix!(int, rowMajor)(3, 5);
    auto b = Matrix!(int, columnMajor)(5, 3);
    assert (a.payload == b.payload);
    a._rawPayload[0][1] = 1;
    b._rawPayload[0][1] = 1;
    assert (a.payload[0][1] == 1);
    assert (b.payload[0][1] == 1);
    a = a * 3;
    b = 3 * b;
    assert (a.payload == b.payload);
    assert (a.row(2).payload == b.column(2).payload);
    assert (!a.byRow.empty);
    assert (!b.byColumn.empty);
    foreach (p; a.byRow.zip(b.byColumn))
        assert (p[0].payload == p[1].payload);
    assert (a.getColumn(4).payload == b.getRow(4).payload);
}

/** Diagonal matrix type.

Scaler multiplication for rows or columns should be translated to the product with Diagonal.
*/
struct Diagonal(T)
    if (is (T == Unqual!T))
{
    /// Initialize Diagonal with payload.
    this (inout T[] payload) inout
    {
        this.payload = payload;
    }
    /// operation as vector
    auto opBinary(string op)(in Diagonal!T rhs) const
    {
        return (Unqual!(typeof (this))(payload.dup)).opOpAssign!op(rhs);
    }
    /// ditto
    auto opOpAssign(string op)(in Diagonal!T rhs)
        if (op == "*" || op == "/" || op == "+" || op == "-")
    {
        mixin ("payload[] "~op~"= rhs.payload[];");
        return this;
    }

    /// scaler multiplication.
    auto opOpAssign(string op)(in T rhs)
        if (op == "*" || op == "/")
    {
        mixin ("payload[] "~op~"= rhs;");
        return this;
    }
    /// ditto
    auto opBinary(string op)(in T rhs) const
        if (op == "*" || op == "/")
    {
        return copy.opOpAssign!op(rhs);
    }
    /// ditto
    auto opBinaryRight(string op)(in T lhs) const
        if (op == "*")
    {
        return opBinary!op(lhs);
    }
    /// ditto
    auto opBinaryRight(string op)(in T lhs) const
        if (op == "/")
    {
        auto payload = this.payload.dup;
        payload[] = lhs / payload[];
        return Unqual!(typeof (this))(payload);
    }

    /// multiplication with vector.
    auto opBinaryRight(string op)(in RowVector!T lhs) const
        if (op == "*")
    {
        auto ret = lhs.copy;
        ret.payload[] *= payload[];
        return ret;
    }
    /// ditto
    auto opBinary(string op)(in ColumnVector!T rhs) const
        if (op == "*")
    {
        auto ret = rhs.copy;
        ret.payload[] *= payload[];
        return ret;
    }

    /// multiplication with matrix
    auto opBinary(string op)(in Matrix!(T, rowMajor) matrix) const
    {
        auto payload = new T[][](matrix.majorLength, matrix.minorLength);
        foreach (i, row; matrix.payload)
            payload[i][] = row[] * this.payload[i];
        return Matrix!(T, rowMajor)(payload);
    }
    /// ditto
    auto opBinary(string op)(in Matrix!(T, columnMajor) matrix) const
    {
        auto payload = new T[][](matrix.majorLength, matrix.minorLength);
        foreach (i, column; matrix.payload)
            payload[i][] = this.payload[] * column[];
        return Matrix!(T, columnMajor)(payload);
    }
    /// ditto
    auto opBinaryRight(string op)(in Matrix!(T, columnMajor) matrix) const
    {
        auto payload = new T[][](matrix.majorLength, matrix.minorLength);
        foreach (i, column; matrix.payload)
            payload[i][] = this.payload[i] * column[];
        return Matrix!(T, columnMajor)(payload);
    }
    /// ditto
    auto opBinaryRight(string op)(in Matrix!(T, rowMajor) matrix) const
    {
        auto payload = new T[][](matrix.majorLength, matrix.minorLength);
        foreach (i, row; matrix.payload)
            payload[i][] = row[] * this.payload[];
        return Matrix!(T, rowMajor)(payload);
    }
package:
    T[] payload;
    auto copy() const
    {
        return Unqual!(typeof (this))(payload.dup);
    }
}

/// ditto
inout(Diagonal!T) diagonal(T)(inout T[] payload)
    if (is (T == Unqual!T))
{
    return inout Diagonal!T(payload);
}
unittest
{
    int[] x;
    const int[] y;
    const(int)[] z;
    static assert (is (typeof (x.diagonal) == Diagonal!int));
    static assert (is (typeof (y.diagonal) == const Diagonal!int));
    static assert (is (typeof (z.diagonal) == const Diagonal!int));
}

/// ditto
auto diagonal(V)(inout V vector)
    if (is (V : RowVector!T, T) ||
        is (V : ColumnVector!T, T))
{
    return vector.payload.diagonal;
}

unittest
{
    assert (([1, 2, 3].row * [1, 2, 4].diagonal).payload == [1, 4, 12]);
    assert (([1, 2, 4].diagonal * [1, 2, 3].column).payload == [1, 4, 12]);
    static assert (!is (typeof ([1].diagonal * [1].row)));
    static assert (!is (typeof ([1].column * [1].diagonal)));
}
unittest
{   import std.string : format;
    auto a = [1, 3, 5].row;
    auto b = [2, 4, 6].row;
    auto c = [1, 2, 4].diagonal;
    assert (a * (b.diagonal + c) == a * b.diagonal + a * c);
    assert (a * (b.diagonal - c) == a * b.diagonal - a * c);
    assert ((a + b) * c == a * c + b * c);
    assert ((a - b) * c == a * c - b * c);
}
unittest
{
    auto a = [1, 3, 5].diagonal;
    auto b = [2, 4, 6].column.diagonal;
    auto c = [1, 2, 4].column;
    assert ((a + b) * c == a * c + b * c);
    assert ((a - b) * c == a * c - b * c);
}
unittest
{
    auto a = [[0, 1, 2], [3, 4, 5]].matrix;
    auto b = [2, 3, 4].diagonal;
    assert ((a * b).payload == [[0, 3, 8], [6, 12, 20]]);
    auto c = [2, 3].diagonal;
    assert ((c * a).payload == [[0, 2, 4], [9, 12, 15]]);

    auto t = [[0, 1, 2], [3, 4, 5]].matrix!columnMajor;
    assert ((b * t).payload == [[0, 3, 8], [6, 12, 20]]);
    assert ((t * c).payload == [[0, 2, 4], [9, 12, 15]]);
}

unittest
{
    auto a = [2, 4, 6, 8].diagonal;
    a = a / 2;
    assert (a.payload == [1, 2, 3, 4]);
    a = 12 / a;
    assert (a.payload == [12, 6, 4, 3]);
    a = 5 * a;
    assert (a.payload == [60, 30, 20, 15]);
}

/// Matrix.opOpAssignRight
void multiplyFromLeft(T, MatrixType)(MatrixType matrix, in Diagonal!T diagonal)
    if (is (MatrixType == Matrix!(T, rowMajor)))
{
    foreach (i, row; matrix.payload)
        row[] *= diagonal.payload[i];
}
/// ditto
void multiplyFromLeft(T, MatrixType)(MatrixType matrix, in Diagonal!T diagonal)
    if (is (MatrixType == Matrix!(T, columnMajor)))
{
    foreach (row; matrix.payload)
        row[] *= diagonal.payload[];
}

unittest
{
    auto a = Matrix!(int, rowMajor)([[1, 2], [3, 4]]);
    auto b = Matrix!(int, columnMajor)([[1, 3], [2, 4]]);
    const c = [5, 6];
    a.multiplyFromLeft(c.diagonal);
    assert (a.payload == [[5, 10], [18, 24]]);
    b.multiplyFromLeft(c.diagonal);
    assert (b.payload == [[5, 18], [10, 24]]);
}

private T[][] minor(T)(in T[][] arr, size_t imin, size_t imax, size_t jmin, size_t jmax)
in
{
    assert (imax <= arr.length);
    if (arr.length)
        assert (jmax <= arr[0].length);
    assert (imin <= imax);
    assert (jmin <= jmax);
}
body
{
    T[][] ret;
    foreach (i; imin..imax)
        ret ~= arr[i][jmin .. jmax].dup;
    return ret;
}
unittest
{
    auto a = [[0, 1, 2], [3, 4, 5], [6, 7, 8]];
    assert (a.minor(0, 3, 0, 3) == a);
    assert (a.minor(0, 2, 0, 2) == [[0, 1], [3, 4]]);
    assert (a.minor(2, 3, 2, 3) == [[8]]);
}

unittest
{
    import std.stdio;
    stderr.writeln("linear.matrix: All green!");
}
