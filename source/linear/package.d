module linear;
public import linear.matrix;
public import linear.vector;
public import linear.saveload;
public import linear.svd;

package import linear.internal;

package import std.range, std.algorithm, std.array, std.experimental.logger;
package import std.string : format;
package import std.traits : Unqual, isFloatingPoint, isIntegral;
