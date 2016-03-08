module linear.internal;

import std.algorithm, std.traits, std.range, std.math, std.array;

package:

bool allFinite(T)(in T[] arr)
    if (is (Unqual!T == T) && isFloatingPoint!T)
{
    foreach (elem; arr)
        if (!elem.isFinite)
            return false;
    return true;
}
bool allFinite(T)(in T[][] arr)
    if (is (Unqual!T == T) && isFloatingPoint!T)
{
    foreach (row; arr)
        if (!row.allFinite)
            return false;
    return true;
}

unittest
{
    assert ([0.0, 1.0].allFinite);
    assert ([[0.0, 1.0], [2.0, 3.0]].allFinite);
    assert (![real.nan].allFinite);
    assert (![real.infinity].allFinite);
}

unittest
{
    import std.stdio;
    stderr.writefln("linear.internal: All green!");
}
