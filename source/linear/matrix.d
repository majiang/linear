module linear.matrix;

import linear;
import std.typecons;
alias RowMajor = Flag!"RowMajor";
enum rowMajor = RowMajor.yes, columnMajor = RowMajor.no;

/** Matrix type. The layout can be chosen from row-major and column-major.
*/
struct Matrix(T, RowMajor rowMajor=rowMajor)
{
	/// Initialize Matrix with given payload.
	this (T[][] payload)
	{
		this.payload = payload;
	}
	/// Zero-clearer.
	static if (T.init != 0)
	void zero()
	{
		foreach (row; payload)
			row[] = 0;
	}
	/// Initialize zero Matrix with specified rows and columns;
	static if (rowMajor)
	this (size_t rows, size_t columns)
	{
		payload = new T[][](rows, columns);
		static if (T.init != 0)
			zero();
	}
	else
	this (size_t rows, size_t columns)
	{
		payload = new T[][](columns, rows);
		static if (T.init != 0)
			zero();
	}
	/// Addition and subtraction.
	auto opOpAssign(string op)(typeof (this) rhs)
		if (op == "+" || op == "-")
	{
		foreach (i, major; rhs.payload)
	mixin ("payload[i][] "~op~"= major[];");
		return this;
	}
	/// ditto
	auto opBinary(string op)(typeof (this) rhs)
		if (op == "+" || op == "-")
	{
		T[][] payload;
		foreach (major; this.payload)
			payload ~= major.dup;
mixin ("return (typeof (this)(payload)) "~op~"= rhs;");
	}
	/// Scaler multiplication.
	auto opOpAssign(string op)(T rhs)
		if (op == "*" || op == "/")
	{
		foreach (major; payload)
	mixin ("major[] "~op~"= rhs;");
		return this;
	}
	/// ditto
	auto opBinary(string op)(T rhs)
		if (op == "*" || op == "/")
	{
		return copy.opOpAssign!op(rhs);
	}
	/// ditto
	auto opBinaryRight(string op)(T lhs)
		if (op == "*")
	{
		return opBinary!op(lhs);
	}
	/// Multiplication with vector.
	auto opBinaryRight(string op)(RowVector!T lhs)
		if (op == "*" && !rowMajor)
	{
		RowVector!T ret;
		foreach (column; this.payload)
			ret.payload ~= lhs * column;
		return ret;
	}
	/// ditto
	auto opBinaryRight(string op)(RowVector!T lhs)
		if (op == "*" && rowMajor)
	{
		auto ret = new T[](this.minorLength()).row;
		static if (T.init != 0)
			ret.payload[] = 0;
		foreach (i, row; this.payload)
			ret.payload[] += lhs.payload[i] * row[];
		return ret;
	}
	/// ditto
	auto opBinary(string op)(ColumnVector!T rhs)
		if (op == "*" && rowMajor)
	{
		ColumnVector!T ret;
		foreach (row; this.payload)
			ret.payload ~= row * rhs;
		return ret;
	}
	/// ditto
	auto opBinary(string op)(ColumnVector!T rhs)
		if (op == "*" && !rowMajor)
	{
		auto ret = new T[](this.minorLength()).column;
		static if (T.init != 0)
			ret.payload[] = 0;
		foreach (i, column; this.payload)
			ret.payload[] += column[] * rhs.payload[i];
		return ret;
	}
	/// Matrix multiplication.
	auto opBinary(string op)(typeof (this) rhs)
		if (op == "*" && rowMajor)
	{
		typeof (this) ret;
		foreach (row; payload)
			ret.payload ~= (RowVector!T(row) * rhs).payload;
		return ret;
	}
	/// ditto
	auto opBinary(string op)(typeof (this) rhs)
		if (op == "*" && !rowMajor)
	{
		typeof (this) ret;
		foreach (column; rhs.payload)
			ret.payload ~= (this * ColumnVector!T(column)).payload;
		return ret;
	}
	auto rows()
	{
		static if (rowMajor)
			return majorLength;
		else
			return minorLength;
	}
	auto columns()
	{
		static if (rowMajor)
			return minorLength;
		else
			return majorLength;
	}
	/// Get a row/column.
	static if (rowMajor)
	auto row(size_t index)
	{
		return payload[index].row;
	}
	else
	auto column(size_t index)
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
	else
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
	auto getColumn(size_t index)
	{
		T[] ret;
		foreach (row; payload)
			ret ~= row[index];
		return ret.column;
	}
	else
	auto getRow(size_t index)
	{
		T[] ret;
		foreach (column; payload)
			ret ~= column[index];
		return ret.row;
	}
	/// deep copy.
	auto copy()
	{
		T[][] ret;
		foreach (major; payload)
			ret ~= major.dup;
		return typeof (this)(ret);
	}
	/// For element-wise special operation.
	auto _rawPayload()
	{
		return payload;
	}
	import std.format;
	void toString(scope void delegate(const(char)[]) sink, FormatSpec!char fmt) const
	{
		sink.formatValue(rowMajor ? payload : payload.transpose, fmt);
	}
package:
	T[][] payload;
private:
	size_t majorLength()
	{
		return payload.length;
	}
	size_t minorLength()
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
	auto a = [[0, 0, 0], [3, 3, 3]].matrix;
	auto b = [[0, 1, 2], [0, 1, 2]].matrix;
	assert ((a + b).payload == [[0, 1, 2], [3, 4, 5]]);
	assert ((a - b).payload == [[0,-1,-2], [3, 2, 1]]);
}
unittest
{
	auto a = [[0, 0, 0], [3, 3, 3]].matrix!columnMajor;
	auto b = [[0, 1, 2], [0, 1, 2]].matrix!columnMajor;
	assert ((a + b).payload == [[0, 1, 2], [3, 4, 5]]);
	assert ((a - b).payload == [[0,-1,-2], [3, 2, 1]]);
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
{
	/// operation as vector
	auto opBinary(string op)(Diagonal!T rhs)
	{
		return (typeof (this)(payload.dup)).opOpAssign!op(rhs);
	}
	/// ditto
	auto opOpAssign(string op)(Diagonal!T rhs)
		if (op == "*" || op == "/" || op == "+" || op == "-")
	{
		mixin ("payload[] "~op~"= rhs.payload[];");
		return this;
	}

	/// scaler multiplication.
	auto opOpAssign(string op)(T rhs)
		if (op == "*" || op == "/")
	{
		mixin ("payload[] "~op~"= rhs;");
		return this;
	}
	/// ditto
	auto opBinary(string op)(T rhs)
		if (op == "*" || op == "/")
	{
		return copy.opOpAssign!op(rhs);
	}
	/// ditto
	auto opBinaryRight(string op)(T lhs)
		if (op == "*")
	{
		return opBinary!op(lhs);
	}
	/// ditto
	auto opBinaryRight(string op)(T lhs)
		if (op == "/")
	{
		auto payload = this.payload.dup;
		payload[] = lhs / payload[];
		return typeof (this)(payload);
	}

	/// multiplication with vector.
	auto opBinaryRight(string op)(RowVector!T lhs)
		if (op == "*")
	{
		auto ret = lhs.copy;
		ret.payload[] *= payload[];
		return ret;
	}
	/// ditto
	auto opBinary(string op)(ColumnVector!T rhs)
		if (op == "*")
	{
		auto ret = rhs.copy;
		ret.payload[] *= payload[];
		return ret;
	}

	/// multiplication with matrix
	auto opBinary(string op)(Matrix!(T, rowMajor) matrix)
	{
		auto payload = new T[][](matrix.majorLength, matrix.minorLength);
		foreach (i, row; matrix.payload)
			payload[i][] = row[] * this.payload[i];
		return Matrix!(T, rowMajor)(payload);
	}
	/// ditto
	auto opBinary(string op)(Matrix!(T, columnMajor) matrix)
	{
		auto payload = new T[][](matrix.majorLength, matrix.minorLength);
		foreach (i, column; matrix.payload)
			payload[i][] = this.payload[] * column[];
		return Matrix!(T, columnMajor)(payload);
	}
	/// ditto
	auto opBinaryRight(string op)(Matrix!(T, columnMajor) matrix)
	{
		auto payload = new T[][](matrix.majorLength, matrix.minorLength);
		foreach (i, column; matrix.payload)
			payload[i][] = this.payload[i] * column[];
		return Matrix!(T, columnMajor)(payload);
	}
	/// ditto
	auto opBinaryRight(string op)(Matrix!(T, rowMajor) matrix)
	{
		auto payload = new T[][](matrix.majorLength, matrix.minorLength);
		foreach (i, row; matrix.payload)
			payload[i][] = row[] * this.payload[];
		return Matrix!(T, rowMajor)(payload);
	}
package:
	T[] payload;
	auto copy()
	{
		return typeof (this)(payload.dup);
	}
}

/// ditto
auto diagonal(T)(T[] payload)
{
	return Diagonal!T(payload);
}
/// ditto
auto diagonal(V)(V vector)
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
{	import std.string : format;
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

unittest
{
	import std.stdio;
	stderr.writeln("linear.matrix: All green!");
}
