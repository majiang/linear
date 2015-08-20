module linear.vector;
import linear;

private mixin template VectorCommon()
{
	auto opIndex()
	{
		struct Result
		{
			auto opOpAssign(string op)(T[] rhs)
			{
				mixin ("payload[] "~op~"= rhs[];");
				return this;
			}
			auto opOpAssign(string op)(Result rhs)
			{
				return opOpAssign!op(rhs.payload);
			}
			auto opBinary(string op)(Result rhs)
			{
				return Result(payload.dup).opOpAssign!op(rhs);
			}
		private:
			T[] payload;
		}
		return Result(payload);
	}
}

struct RowVector(T)
{
	mixin VectorCommon;
	auto opOpAssign(string op)(RowVector!T rhs)
		if (op == "+" || op == "-")
	{
		pragma (msg, typeof (this), ".opOpAssign!", op);
		mixin ("payload[] "~op~"= rhs.payload[];");
		return this;
	}
	auto opBinary(string op)(RowVector!T rhs)
		if (op == "+" || op == "-")
	{
		pragma (msg, typeof (this), ".opBinary!", op);
		return this.copy.opOpAssign!op(rhs);
	}
	auto opBinary(string op)(ColumnVector!T rhs)
		if (op == "*")
	{
		return this.product(rhs);
	}
package:
	T[] payload;
	auto copy()
	{
		return typeof (this)(payload.dup);
	}
}
auto row(T)(T[] payload)
{
	return RowVector!T(payload);
}
struct ColumnVector(T)
{
	mixin VectorCommon;
	auto opOpAssign(string op)(typeof (this) rhs)
		if (op == "+" || op == "-")
	{
		pragma (msg, typeof (this), ".opOpAssign!", op);
		mixin ("payload[] "~op~"= rhs.payload[];");
		return this;
	}
	auto opBinary(string op)(typeof (this) rhs)
	{
		pragma (msg, typeof (this), ".opBinary!", op);
		return this.copy.opOpAssign!op(rhs);
	}
package:
	T[] payload;
	auto copy()
	{
		return typeof (this)(payload.dup);
	}
}
auto column(T)(T[] payload)
{
	return ColumnVector!T(payload);
}

auto product(T)(RowVector!T x, ColumnVector!T y)
{
	return x.payload.zip(y.payload)
		.map!(a => a[0] * a[1])
		.reduce!((a, b) => a + b);
}

unittest
{
	auto a = [0, 1, 2].row;
	auto b = [3, 4, 5].column;
	assert (a * b == 14);
	assert ((a[] * [3, 4, 5].row[]).payload == [0, 4, 10]);
}

unittest
{
	import std.stdio;
	stderr.writeln("linear.vector: All green!");
}
