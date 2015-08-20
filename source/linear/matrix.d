module linear.matrix;

import linear;
import std.typecons;
alias RowMajor = Flag!"RowMajor";
enum rowMajor = RowMajor.yes, columnMajor = RowMajor.no;

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


struct Diagonal(T)
{
	auto opBinary(string op)(Diagonal!T rhs)
	{
		return (typeof (this)(payload.dup)).opOpAssign!op(rhs);
	}
	auto opOpAssign(string op)(Diagonal!T rhs)
		if (op == "*" || op == "/" || op == "+" || op == "-")
	{
		mixin ("payload[] "~op~"= rhs.payload[];");
		return this;
	}
	auto opBinaryRight(string op)(RowVector!T lhs)
		if (op == "*")
	{
		auto ret = lhs.copy;
		ret.payload[] *= payload[];
		return ret;
	}
	auto opBinary(string op)(ColumnVector!T rhs)
		if (op == "*")
	{
		auto ret = rhs.copy;
		ret.payload[] *= payload[];
		return ret;
	}
private:
	T[] payload;
package:
	auto copy()
	{
		return typeof (this)(payload.dup);
	}
}

auto diagonal(T)(T[] payload)
{
	return Diagonal!T(payload);
}
auto diagonal(T)(RowVector!T vector)
{
	return vector.payload.diagonal;
}
auto diagonal(T)(ColumnVector!T vector)
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
{
	auto a = [1, 3, 5].row;
	auto b = [2, 4, 6].row;
	auto c = [1, 2, 4].diagonal;
//	assert (a * (b.diagonal + c) == a * b.diagonal + a * c);
//	assert (a * (b.diagonal - c) == a * b.diagonal - a * c);
//	assert ((a + b) * c == a * b + a * c);
//	assert ((a - b) * c == a * b - a * c);
}
unittest
{
//	auto a = [1, 3, 5].diagonal;
//	auto b = [2, 4, 6].column.diagonal;
//	auto c = [1, 2, 4].column;
//	assert ((a + b) * c == a * c + b * c);
//	assert ((a - b) * c == a * c - b * c);
}

unittest
{
	import std.stdio;
	stderr.writeln("linear.matrix: All green!");
}
