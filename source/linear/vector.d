module linear.vector;
import linear;
import std.traits : Unqual;

/// Row vector.
struct RowVector(T)
	if (is (T == Unqual!T))
{
	/// Zero vector with given length.
	this (in size_t length)
	{
		this.payload = new T[length];
		static if (T.init != 0)
			this.payload[] = 0;
	}
	/// Initialize RowVector with payload.
	this (T[] payload)
	{
		this.payload = payload;
	}
	/// ditto
	this (in T[] payload) const
	{
		this.payload = payload;
	}
	/// Componentwise product.
	auto opIndexOpAssign(string op)(in typeof (this) rhs)
		if (op == "*")
	{
		mixin ("payload[] "~op~"= rhs.payload[];");
		return this;
	}
	/// ditto
	auto opBinary(string op)(in typeof (this) rhs) const
		if (op == "*")
	{
		return this.copy.opIndexOpAssign!op(rhs);
	}

	/// Vector addition and subtraction.
	auto opOpAssign(string op)(in typeof (this) rhs)
		if (op == "+" || op == "-")
	{
		mixin ("payload[] "~op~"= rhs.payload[];");
		return this;
	}
	/// ditto
	auto opBinary(string op)(in typeof (this) rhs) const
		if (op == "+" || op == "-")
	{
		return this.copy.opOpAssign!op(rhs);
	}

	/// Scalar multiplication and division.
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
		return this.copy.opOpAssign!op(rhs);
	}
	/// ditto
	auto opBinaryRight(string op)(in T lhs) const
		if (op == "*")
	{
		return this.opBinary!op(lhs);
	}

	/// inner product.
	auto opBinary(string op)(in ColumnVector!T rhs) const
		if (op == "*")
	{
		pragma (msg, "Row!T * Column!T");
		return this * rhs.payload;
	}
	/// ditto
	auto opBinary(string op)(in T[] rhs) const
		if (op == "*")
	{
		return this.payload.product(rhs);
	}
	/// Length.
	auto length() const
	{
		return payload.length;
	}
	/// For special element-wise operation.
	auto opIndex() const
	{
		return payload;
	}
	/// ditto
	auto opIndex()
	{
		return payload;
	}
	/// Deep copy.
	auto copy() const
	{
		return Unqual!(typeof (this))(payload.dup);
	}
	import std.format;
	void toString(scope void delegate(const(char)[]) sink, FormatSpec!char fmt) const
	{
		sink.formatValue(payload, fmt);
	}
package:
	T[] payload;
}
/// ditto
auto row(T)(in T[] payload)
{
	return const RowVector!T(payload);
}
/// ditto
auto row(T)(T[] payload)
	if (is (T == Unqual!T))
{
	return RowVector!T(payload);
}
unittest
{
	int[] x;
	const int[] y;
	const(int)[] z;
	static assert (is (typeof (x.row) == RowVector!int));
	static assert (is (typeof (y.row) == const RowVector!int));
	static assert (is (typeof (z.row) == const RowVector!int));
}

/// Column vector.
struct ColumnVector(T)
	if (is (T == Unqual!T))
{
	/// Zero vector with given length.
	this (in size_t length)
	{
		this.payload = new T[length];
		static if (T.init != 0)
			this.payload[] = 0;
	}
	/// Initialize ColumnVector with payload.
	this (T[] payload)
	{
		this.payload = payload;
	}
	/// Initialize ColumnVector with payload.
	this (in T[] payload) const
	{
		this.payload = payload;
	}
	/// Componentwise product.
	auto opIndexOpAssign(string op)(in typeof (this) rhs)
		if (op == "*")
	{
		mixin ("payload[] "~op~"= rhs.payload[];");
		return this;
	}
	/// ditto
	auto opBinary(string op)(in typeof (this) rhs) const
		if (op == "*")
	{
		return this.copy.opIndexOpAssign!op(rhs);
	}

	/// Vector addition and subtraction.
	auto opOpAssign(string op)(in typeof (this) rhs)
		if (op == "+" || op == "-")
	{
		mixin ("payload[] "~op~"= rhs.payload[];");
		return this;
	}
	/// ditto
	auto opBinary(string op)(in typeof (this) rhs) const
		if (op == "+" || op == "-")
	{
		return this.copy.opOpAssign!op(rhs);
	}

	/// Scalar multiplication and division.
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
		return this.copy.opOpAssign!op(rhs);
	}
	/// ditto
	auto opBinaryRight(string op)(in T lhs) const
		if (op == "*")
	{
		return this.opBinary!op(lhs);
	}

	/// inner product.
	auto opBinaryRight(string op)(in T[] lhs) const
		if (op == "*")
	{
		return lhs.product(this.payload);
	}
	/// Length.
	auto length() const
	{
		return payload.length;
	}
	/// For special element-wise operation.
	auto opIndex() const
	{
		return payload;
	}
	/// ditto
	auto opIndex()
	{
		return payload;
	}
	/// Deep copy.
	auto copy() const
	{
		return Unqual!(typeof (this))(payload.dup);
	}
	import std.format;
	void toString(scope void delegate(const(char)[]) sink, FormatSpec!char fmt) const
	{
		sink.formatValue(payload, fmt);
	}
package:
	T[] payload;
}
/// ditto
auto column(T)(in T[] payload)
{
	return const ColumnVector!T(payload);
}
/// ditto
auto column(T)(T[] payload)
	if (is (T == Unqual!T))
{
	return ColumnVector!T(payload);
}
unittest
{
	int[] x;
	const int[] y;
	const(int)[] z;
	static assert (is (typeof (x.column) == ColumnVector!int));
	static assert (is (typeof (y.column) == const ColumnVector!int));
	static assert (is (typeof (z.column) == const ColumnVector!int));
}

/// fill by 1s.
auto ones(T)(size_t length)
{
	return T(1).repeat.take(length).array;
}
/// transposition of vector.
auto transpose(T)(in RowVector!T row)
{
	return row.payload.column;
}
/// ditto
auto transpose(T)(RowVector!T row)
{
	return row.payload.column;
}
/// ditto
auto transpose(T)(in ColumnVector!T column)
{
	return column.payload.row;
}
/// ditto
auto transpose(T)(ColumnVector!T column)
{
	return column.payload.row;
}

/// inner product
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
	auto a = [0, 2, 4];
	auto r = a.row;
	auto c = a.column;
	assert (r[] is c[]);
	assert (r.length == c.length);
}
unittest
{
	auto a = [1, 2, 3].row;
	auto b = a.length.ones!int.column;
	assert (b.transpose * a.transpose == a * b);
	assert ((cast(const)b).transpose * (cast(const)a).transpose == a * b);
}

unittest
{
	import std.stdio;
	stderr.writeln("linear.vector: All green!");
}
