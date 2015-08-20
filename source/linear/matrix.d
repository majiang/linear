module linear.matrix;

import linear;
import std.typecons;
alias RowMajor = Flag!"RowMajor";
enum rowMajor = RowMajor.yes, columnMajor = RowMajor.no;

/** Matrix type. The layout can be chosen from row-major and column-major.
*/
struct Matrix(T, RowMajor rowMajor=rowMajor)
{
	auto opOpAssign(string op)(typeof (this) rhs)
		if (op == "+" || op == "-")
	{
		foreach (i, major; rhs.payload)
	mixin ("payload[i][] "~op~"= major[];");
		return this;
	}
	auto opBinary(string op)(typeof (this) rhs)
		if (op == "+" || op == "-")
	{
		T[][] payload;
		foreach (major; this.payload)
			payload ~= major.dup;
mixin ("return (typeof (this)(payload)) "~op~"= rhs;");
	}
	auto opBinaryRight(string op)(RowVector!T lhs)
		if (op == "*" && !rowMajor)
	{
		RowVector!T ret;
		foreach (column; this.payload)
			ret.payload ~= lhs * column;
		return ret;
	}
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
	auto opBinary(string op)(ColumnVector!T rhs)
		if (op == "*" && rowMajor)
	{
		ColumnVector!T ret;
		foreach (row; this.payload)
			ret.payload ~= row * rhs;
		return ret;
	}
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
	auto opBinary(string op)(typeof (this) rhs)
		if (op == "*" && rowMajor)
	{
		typeof (this) ret;
		foreach (row; payload)
			ret.payload ~= (RowVector!T(row) * rhs).payload;
		return ret;
	}
	auto opBinary(string op)(typeof (this) rhs)
		if (op == "*" && !rowMajor)
	{
		typeof (this) ret;
		foreach (column; rhs.payload)
			ret.payload ~= (this * ColumnVector!T(column)).payload;
		return ret;
	}
private:
	T[][] payload;
	size_t majorLength()
	{
		return payload.length;
	}
	size_t minorLength()
	{
		return majorLength ? payload.front.length : 0;
	}
}

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

/* Diagonal matrix type.

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

auto diagonal(T)(T[] payload)
{
	return Diagonal!T(payload);
}
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
	import std.stdio;
	stderr.writeln("linear.matrix: All green!");
}
