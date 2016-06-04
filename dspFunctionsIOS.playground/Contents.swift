//: Playground - noun: a place where people can play
//This IOS playground is intended to be a reference for DSP functions in the Accelerate framework. A struct is used for vectors and matrices definition.


//5/31/2016
import UIKit
import Accelerate
import XCPlayground

var str = "Hello, playground"



public struct Vector {
    let rows :Int
    var grid=[Double]()
    init(rows: Int) {
        self.rows = rows
        grid = Array(count: rows, repeatedValue: 0.0)
    }
    func indexIsValid(row: Int) -> Bool {
        return row >= 0 && row < rows
    }
    subscript(row: Int) -> Double {
        get {
            assert(indexIsValid(row), "Index out of range")
            return grid[row]
        }
        set {
            assert(indexIsValid(row), "Index out of range")
            grid[row] = newValue
        }
    }
}


public struct Matrix {
    let rows: Int, columns: Int
    var grid: [Double]
    
    init(rows: Int, columns: Int) {
        self.rows = rows
        self.columns = columns
        grid = Array(count: rows * columns, repeatedValue: 0.0)
    }
    
    func indexIsValidForRow(row: Int, column: Int) -> Bool {
        return row >= 0 && row < rows && column >= 0 && column < columns
    }
    
    subscript(row: Int, column: Int) -> Double {
        get {
            assert(indexIsValidForRow(row,column: column), "Index out of range")
            return grid[(row * columns) + column]
        }
        set {
            assert(indexIsValidForRow(row, column: column), "Index out of range")
            grid[(row * columns) + column] = newValue
        }
    }
}


// A little playground helper function
func plotArrayInPlayground<T>(arrayToPlot:Array<T>, title:String) {
    for currentValue in arrayToPlot {
        XCPCaptureValue(title, value: currentValue)
    }
}

func abs (a: [Double]) -> [Double] {
    //Vector absolute values; double precision.
    var result = [Double](count:a.count, repeatedValue:0.0)
    vDSP_vabsD(a, 1,&result, 1, UInt(a.count))
    return result
}

func ramp (initial_value: Double, increment: Double, numberOfSteps:Int) -> [Double] {
    //Build ramped vector; double precision.
    var aa = initial_value
    var bb=increment
    
    var result = [Double](count:numberOfSteps, repeatedValue:0.0)
    vDSP_vrampD(&aa, &bb, &result, 1, UInt(result.count))
    return result
}

func neg (a: [Double]) -> [Double] {
    //Computes the squared values of vector A and leaves the result in vector C; double precision.
    var result = [Double](count:a.count, repeatedValue:0.0)
    vDSP_vnegD(a, 1,&result, 1, UInt(a.count))
    return result
}

func sq (a: [Double]) -> [Double] {
    //Computes the squared values of vector A and leaves the result in vector C; double precision.
    var result = [Double](count:a.count, repeatedValue:0.0)
    vDSP_vsqD(a, 1,&result, 1, UInt(a.count))
    return result
}

func normalize (a: [Double]) -> (out:[Double], mean1:Double, std:Double) {
    //Compute mean and standard deviation and then calculate new elements to have a zero mean and a unit standard deviation. Double precision.
    var mean = 0.0
    var standardDeviation = 0.0
    
    var result = [Double](count:a.count, repeatedValue:0.0)
    vDSP_normalizeD(a, 1,&result, 1, &mean, &standardDeviation, UInt(a.count))
    return (out:result, mean1:mean, std:standardDeviation)
}

func polar (a: [Double]) -> [Double] {
    //Converts rectangular coordinates to polar coordinates. Cartesian (x,y) pairs are read from vector A. Polar (rho, theta) pairs, where rho is the radius(magnitude) and theta is the angle in radians in the range [-pi, pi] are written to vector C. N specifies the number of coordinate pairs in A and C. Coordinate pairs are adjacent elements in the array, regardless of stride; stride is the distance from one coordinate pair to the next.
    var result = [Double](count:a.count, repeatedValue:0.0)
    vDSP_polarD(a,2,&result,2, UInt(a.count/2))
    return result
}

func rect (a: [Double]) -> [Double] {
    //Converts polar coordinates to rectangular coordinates. Polar (rho, theta) pairs, where rho is the radius and theta is the angle in the range [-pi, pi] are read from vector A. Cartesian (x,y) pairs are written to vector C. N specifies the number of coordinate pairs in A and C.
    //Coordinate pairs are adjacent elements in the array, regardless of stride; stride is the distance from one coordinate pair to the next.
    var result = [Double](count:a.count, repeatedValue:0.0)
    vDSP_rectD(a, 2,&result, 2, UInt(a.count/2))
    return result
}

func add (a: [Double], b: [Double]) -> [Double] {
    //This function adds the first N elements of A to B and leaves the result in C:

    assert(a.count == b.count, "Expected arrays of the same length, instead got arrays of two different lengths")
    
    var result = [Double](count:a.count, repeatedValue:0.0)
    vDSP_vaddD(a, 1, b, 1, &result, 1, UInt(a.count))
    return result
}

func sum (a: [Double]) -> Double {
    //Takes the sum of the input array.
    var result = 0.0
    vDSP_sveD(a, 1, &result, UInt(a.count))
    return result
}

func sub (a: [Double], b: [Double]) -> [Double] {
    //This function subtracts the first N elements of A - B and leaves the result in C. This function subtracts the first N elements of B from A and leaves the result in C.
    
    assert(a.count == b.count, "Expected arrays of the same length, instead got arrays of two different lengths")
    
    var result = [Double](count:a.count, repeatedValue:0.0)
    vDSP_vsubD(a, 1, b, 1, &result, 1, UInt(a.count))
    return result
}

func mul (a: [Double], b: [Double]) -> [Double] {
    //This function multiplies the first N elements of A by corresponding elements of B, leaving the result in C.
    assert(a.count == b.count, "Expected arrays of the same length, instead got arrays of two different lengths")
    var result = [Double](count:a.count, repeatedValue:0.0)
    vDSP_vmulD(a, 1, b, 1, &result, 1, UInt(a.count))
    return result
}

func div (a: [Double], b: [Double]) -> [Double] {
    //This function divides the first N elements of A by corresponding elements of B, leaving the result in C.
    assert(a.count == b.count, "Expected arrays of the same length, instead got arrays of two different lengths")
    var result = [Double](count:a.count, repeatedValue:0.0)
    vDSP_vdivD(a, 1, b, 1, &result, 1, UInt(a.count))
    return result
}

func dotProduct (a: [Double], b: [Double]) -> Double {
    //Computes the dot or scalar product of vectors A and B and leaves the result in scalar *C; double precision.
    assert(a.count == b.count, "Expected arrays of the same length, instead got arrays of two different lengths")
    var result = 0.0
    vDSP_dotprD(a, 1, b, 1, &result, UInt(a.count))
    return result
}


func transpose (a: Matrix ) -> [Double] {
    //Creates a transposed matrix C from a source matrix A; double precision.
    let rows=a.rows
    let columns=a.columns
    var a1=[Double]()
    var count=0
    a1 = [Double](count:rows*columns, repeatedValue:0.0)
    for i in 0..<columns{
        for j in 0..<rows{
             print("\(count)")
            a1[count]=a[i,j]
            count += 1
    } }
    var result = [Double](count:rows*columns, repeatedValue:0.0)
    vDSP_mtransD(a1, 1, &result, 1, UInt(rows),UInt(columns))
    return result
}

func mmul (a: Matrix, b: Matrix) -> [Double]{
    //This function multiplies an M-by-P matrix A by a P-by-N matrix B and stores the results in an M-by-N matrix C.
    
    let rows=a.rows
    let columns=a.columns
    //let rows1=b.rows
    let columns1=b.columns
    
    var a1=[Double]()
    var b1=[Double]()
    var count=0
    
    a1 = [Double](count:rows*columns, repeatedValue:0.0)
    for i in 0..<columns{
        for j in 0..<rows{
            print("\(count)")
            a1[count]=a[i,j]
            count += 1
        }
    }
    print("\(a1)")
    
    b1 = [Double](count:rows*columns, repeatedValue:0.0)
    
    count=0
    for i in 0..<columns{
        for j in 0..<rows{
            print("\(count)")
            b1[count]=b[i,j]
            count += 1
        }
    }
    print("\(b1)")
    
    var result = [Double](count:rows*columns, repeatedValue:0.0)
    //func vDSP_mmulD(_ __A: UnsafePointer<Double>, _ __IA: vDSP_Stride, _ __B: UnsafePointer<Double>, _ __IB: vDSP_Stride, _ __C: UnsafeMutablePointer<Double>, _ __IC: vDSP_Stride, _ __M: vDSP_Length, _ __N: vDSP_Length, _ __P: vDSP_Length)
    
    vDSP_mmulD(a1, 1, b1, 1, &result,1, UInt(rows),UInt(columns1),UInt(rows))
    return result
}


func distance (a: [Double], b: [Double]) -> [Double] {
//For the first N elements of A and B, this function takes the square root of the sum of the squares of corresponding elements, leaving the results in C. Vector distance; double precision.
    assert(a.count == b.count, "Expected arrays of the same length, instead got arrays of two different lengths")
    var result = [Double](count:a.count, repeatedValue:0.0)
    vDSP_vdistD(a, 1, b, 1, &result, 1, UInt(a.count))
    return result
}

func mvmul (a: Matrix, b: Vector) -> [Double]{
    //This function multiplies an M-by-P matrix A by a M-by-1 column vector B and stores the results in an M-by-1 vector C.
    
    let rows=a.rows
    let columns=a.columns
    let rows1=b.rows
    
    var a1=[Double]()
    var b1=[Double]()
    var count=0
    
    a1 = [Double](count:rows*columns, repeatedValue:0.0)
    for i in 0..<columns{
        for j in 0..<rows{
            print("\(count)")
            a1[count]=a[i,j]
            count += 1
        }
    }
    print("\(a1)")
    
    b1 = [Double](count:rows1, repeatedValue:0.0)
    
    count=0
    for i in 0..<rows1{
            print("\(count)")
            b1[count]=b[i]
            count += 1
    }
    print("\(b1)")
    
    var result = [Double](count:rows1, repeatedValue:0.0)
    
    for i in 0..<columns{
        for j in 0..<rows1{
        result[i] += a[i,j]*b[j]
        }
    }
    
    return result
}

func gauss(a:Matrix)->[Double]{
    //gauss elimination to solve linear equation
    let rows=a.rows
    let columns=a.columns
    var A=Matrix(rows: rows, columns: columns)
    A=a
    var n=rows
    
    var c: Double
    var x: [Double]
    x = Array(count: rows, repeatedValue: 0.0)
    
    // loop for the generation of upper triangular matrix
    /*
     
     for(j=1; j<=n; j++) /* loop for the generation of upper triangular matrix*/
     {
     for(i=1; i<=n; i++)
     {
     if(i>j)
     {
     c=A[i][j]/A[j][j];
     for(k=1; k<=n+1; k++)
     {
     A[i][k]=A[i][k]-c*A[j][k];
     }
     }
     }
     }
     
     */
    // loop for the generation of upper triangular matrix
    //print("\(A)")
    for j in 0..<n{
        for i in 0..<n{
            if i>j{
                c=A[i,j]/A[j,j]
                for k in 0..<n+1{
                    A[i,k] -= c*A[j,k]
                    //print("A[\(i),\(k)]=\(A[i,k])")
                }
            }
        }
    }
    
    //print("\(A)")
    n -= 1
    x[n]=A[n,n+1]/A[n,n]
    //print("x[\(n)]=\(x[n])")
    
    var sum=0.0
    /* this loop is for backward substitution
     for(i=n-1; i>=1; i--)
     {
     sum=0;
     for(j=i+1; j<=n; j++)
     {
     sum=sum+A[i][j]*x[j];
     }
     x[i]=(A[i][n+1]-sum)/A[i][i];
     }
     */
    var i :Int
    
    for(i=n-1; i>=0; i -= 1){
        sum=0.0
        
        for j in i+1..<rows{
            sum += A[i,j]*x[j]
            //print("sum[\(j)]=\(sum)")
            
        }
        x[i]=(A[i,n+1]-sum)/A[i,i]
        //print("x[\(n)]=\(x[n])")
        
    }
    return x
    
}

func max (a: Vector) -> Double {
    //Vector maximum magnitude; double precision.
    var result = 0.0
    var aa=[Double](count:a.rows, repeatedValue:0.0)
    for i in 0..<a.rows{
        aa[i]=a[i]
    }
    vDSP_maxmgvD(aa,1,&result,UInt(a.rows))
    return result
}

public func max1 (a: [Double]) -> Double {
    //Vector maximum magnitude; double precision.
    var result = 0.0
    vDSP_maxmgvD(a,1,&result,UInt(a.count))
    return result
}


var m1 = Matrix(rows: 3, columns: 3)
m1[0, 0] = 1.0
m1[0, 1] = 2.0
m1[0, 2] = 3.0
m1[1, 0] = 3.0
m1[1, 1] = 2.0
m1[1, 2] = 1.0
m1[2, 0] = 2.0
m1[2, 1] = 1.0
m1[2, 2] = 3.0


var solution1=Vector(rows: 3)

solution1[0] = -7.809086
solution1[1] = -8.690904
solution1[2] = 7.418178





