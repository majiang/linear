module linear.saveload;
import linear;
import std.conv, std.string, std.traits;

string toSerialString(Vector)(in Vector v)
	if (is (Vector == RowVector!T, T) ||
		is (Vector == ColumnVector!T, T))
{
	static if (is (Vector == RowVector!T, T) || is (Vector == ColumnVector!T, T))
	{
		static if (isIntegral!T)
			return "[%(%d,%)]".format(v.payload);
		else static if (isFloatingPoint!T)
			return "[%(%a,%)]".format(v.payload);
		else static assert (false, "unknown type");
	}
	else static assert (false, "contradiction");
}

RowVector!T toRowVector(T)(in string s)
{
	return s.to!(T[]).row;
}
ColumnVector!T toColumnVector(T)(in string s)
{
	return s.to!(T[]).column;
}

unittest
{
	auto a = [0, 1, 2].row;
	const b = [3, 4, 5].column;
	immutable c = [real(0), real.infinity, real.epsilon].row;
	assert (a.toSerialString.toRowVector!int == a);
	assert (b.toSerialString.toColumnVector!int == b);
	assert (c.toSerialString.toRowVector!real == c);
}

string toSerialString(T, RowMajor rowMajor)(in Matrix!(T, rowMajor) m)
body
{
	static if (isIntegral!T)
		return "[%([%(%d,%)]%|,%)]".format(m.payload);
	else static if (isFloatingPoint!T)
		return "[%([%(%a,%)]%|,%)]".format(m.payload);
	else static assert (false, "unknown type");
}

auto toMatrix(T, RowMajor rowMajor=rowMajor)(in string s)
{
	return s.to!(T[][]).matrix!rowMajor;
}

unittest
{
	auto a = [[0, 1, 2], [3, 4, 5]].matrix;
	assert (a.toSerialString.toMatrix!int == a);
	immutable c = [[real(0), 1], [real.infinity, real.epsilon]].matrix!columnMajor;
	assert (c.toSerialString.toMatrix!(real, columnMajor) == c);
}

unittest
{
	import std.stdio;
	stderr.writeln("linear.saveload: All green!");
}
