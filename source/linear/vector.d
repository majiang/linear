module linear.vector;
import linear;

struct RowVector(T)
{
	/// Componentwise product.
	auto opIndexOpAssign(string op)(typeof (this) rhs)
		if (op == "*")
	{
		mixin ("payload[] "~op~"= rhs.payload[];");
		return this;
	}
	/// ditto
	auto opBinary(string op)(typeof (this) rhs)
		if (op == "*")
	{
		return this.copy.opIndexOpAssign!op(rhs);
	}

	/// Vector addition and subtraction.
	auto opOpAssign(string op)(typeof (this) rhs)
		if (op == "+" || op == "-")
	{
		mixin ("payload[] "~op~"= rhs.payload[];");
		return this;
	}
	/// ditto
	auto opBinary(string op)(typeof (this) rhs)
		if (op == "+" || op == "-")
	{
		return this.copy.opOpAssign!op(rhs);
	}

	/// Scalar multiplication and division.
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
		return this.copy.opOpAssign!op(rhs);
	}
	/// ditto
	auto opBinaryRight(string op)(T lhs)
		if (op == "*")
	{
		return this.opBinary!op(lhs);
	}

	/// inner product.
	auto opBinary(string op)(ColumnVector!T rhs)
		if (op == "*")
	{
		pragma (msg, "Row!T * Column!T");
		return this * rhs.payload;
	}
	/// ditto
	auto opBinary(string op)(T[] rhs)
		if (op == "*")
	{
		return this.payload.product(rhs);
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
	/// Componentwise product.
	auto opIndexOpAssign(string op)(typeof (this) rhs)
		if (op == "*")
	{
		mixin ("payload[] "~op~"= rhs.payload[];");
		return this;
	}
	/// ditto
	auto opBinary(string op)(typeof (this) rhs)
		if (op == "*")
	{
		return this.copy.opIndexOpAssign!op(rhs);
	}

	/// Vector addition and subtraction.
	auto opOpAssign(string op)(typeof (this) rhs)
		if (op == "+" || op == "-")
	{
		mixin ("payload[] "~op~"= rhs.payload[];");
		return this;
	}
	/// ditto
	auto opBinary(string op)(typeof (this) rhs)
		if (op == "+" || op == "-")
	{
		return this.copy.opOpAssign!op(rhs);
	}

	/// Scalar multiplication and division.
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
		return this.copy.opOpAssign!op(rhs);
	}
	/// ditto
	auto opBinaryRight(string op)(T lhs)
		if (op == "*")
	{
		return this.opBinary!op(lhs);
	}

	/// inner product.
	auto opBinaryRight(string op)(T[] lhs)
		if (op == "*")
	{
		return lhs.product(this.payload);
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

auto product(T)(T[] x, T[] y)
{
	return x.zip(y)
		.map!(a => a[0] * a[1])
		.reduce!((a, b) => a + b);
}

unittest
{
	auto a = [0, 1, 2].row;
	auto b = [3, 4, 5].column;
	assert (a * b == 14);
	assert ((a * [3, 4, 5].row).payload == [0, 4, 10]);
	assert (([0, 1, 2].column * b).payload == [0, 4, 10]);
	assert ((a - 2 * a).payload == [0, -1, -2]);
	assert ((b + 2 * b).payload == [9, 12, 15]);
}

unittest
{
	import std.stdio;
	stderr.writeln("linear.vector: All green!");
}
