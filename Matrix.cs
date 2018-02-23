using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MatrixOp
{
    /// <summary>
    /// This class describes matrix and method's to work with it
    /// </summary>
    class Matrix
    {
        private int n;
        private int m;
        private double[,] matrix;
        Random rand = new Random();


        /// <summary>
        /// creates nxm matrix with random values
        /// </summary>
        /// <param name="n">The amount of raws</param>
        /// <param name="m">The amount of columns</param>
        public Matrix(int n, int m)
        {
            this.n = n;
            this.m = m;
            this.matrix = new double[this.n, this.m];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    this.matrix[i, j] = (double)rand.Next(0, 100);
                }
            }
        }

        /// <summary>
        /// creates nxm matrix with initial values
        /// </summary>
        /// <param name="n">The amount of raws</param>
        /// <param name="m">The amount of columns</param>
        /// <param name="values">The initial values</param>
        public Matrix(int n, int m, double[,] values)
        {
            if (n <= 0 || m <= 0)
                throw new System.ArgumentException("rows and columns must be natural numbers");
            if (n == values.GetLength(0) && m == values.GetLength(1))
            {
                this.n = n;
                this.m = m;
                this.matrix = new double[this.n, this.m];
                for (int i = 0; i < this.n; i++)
                {
                    for (int j = 0; j < this.m; j++)
                    {
                        this.matrix[i, j] = values[i, j];
                    }
                }
            }
            else
            {
                throw new System.ArgumentException("the count of values is incorrect");
            }
        }

        /// <summary>
        /// calculates the addition of two matrix instanses
        /// </summary>
        /// <param name="op1">the left hand side parameter of the operator</param>
        /// <param name="op2">the right hand side parameter of the operator</param>
        /// <returns>the addition of two matrix instanses</returns>
        public static Matrix operator +(Matrix op1, Matrix op2)
        {
            Matrix sum = new Matrix(op1.n, op1.m);
            for (int i = 0; i < sum.n; i++)
            {
                for (int j = 0; j < sum.m; j++)
                {
                    sum.matrix[i, j] = op1.matrix[i, j] + op2.matrix[i, j];
                }
            }
            return sum;
        }

        /// <summary>
        /// calculates the multiplication of two matrix instanses
        /// </summary>
        /// <param name="op1">the left hand side parameter of the operator</param>
        /// <param name="op2">the right hand side parameter of the operator</param>
        /// <returns>the multiplication of two matrix instanses</returns>
        public static Matrix operator *(Matrix op1, Matrix op2)
        {
            if (op1.matrix.GetLength(1) == op2.matrix.GetLength(0))
            {
                Matrix prod = new Matrix(op1.n, op2.m);
                for (int i = 0; i < prod.n; i++)
                {
                    for (int j = 0; j < prod.m; j++)
                    {
                        prod.matrix[i, j] = 0;
                        for (int k = 0; k < op1.m; k++)
                            prod.matrix[i, j] += op1.matrix[i, k] * op2.matrix[k, j];
                    }
                }
                return prod;
            }
            else
            {
                Console.WriteLine("The number of columns in first matrix must be equal to the number of rows in second");
                return null;
            }
        }


       
        /// <summary>
        /// calculates the determinant of the matrix
        /// </summary>
        /// <returns>returns the determinant of the given matrix if it is defined</returns>
        public double Determinant()
        {
            if (this.NullRowCol())
                return 0;

            if (this.n != this.m)
            {
                Console.WriteLine("Matrix determinant is defined for only square matrices, this matrix isn't square.");
                return 0;
            }
            Matrix copyMatrix = new Matrix(this.n, this.m, this.matrix);
            int i, j, k;

            // counts the number of swaps to determine the sign of det
            int count = 0;
            double det = 0;
            for (i = 0; i < copyMatrix.n; i++)
            {
                for (j = i + 1; j < copyMatrix.n; j++)
                {
                    if (copyMatrix.matrix[i, i] == 0)
                    {
                        for (k = i + 1; k < copyMatrix.n; k++)
                        {
                            if (copyMatrix.matrix[k, i] != 0)
                            {
                                double temp;
                                for (int t = 0; t < this.m; t++)
                                {
                                    temp = this.matrix[i, t];
                                    this.matrix[i, t] = this.matrix[k, t];
                                    this.matrix[k, t] = temp;
                                }
                                count++;
                                break;
                            }
                        }
                    }

                    det = copyMatrix.matrix[j, i] / copyMatrix.matrix[i, i];
                    for (k = i; k < copyMatrix.n; k++)
                        copyMatrix.matrix[j, k] = copyMatrix.matrix[j, k] - det * copyMatrix.matrix[i, k];
                }
            }
            det = 1;
            for (i = 0; i < copyMatrix.n; i++)
                det = det * copyMatrix.matrix[i, i];

            //anytime when 2 rows are swapped in matrix determain chanes it's sign so when the count of swaps is odd det is negative

            if (count % 2 == 1)
                det = -det;
            return det;
        }

        /// <summary>
        /// Checks if there is any null(0) row or column in matrix
        /// it means that the determinant of that matrix is 0
        /// </summary>
        /// <returns>True - if there is any null(0) row or column in matrix. False - otherwise</returns>
        public bool NullRowCol()
        {

            int i, j;
            for (i = 0; i < this.n; i++)
            {
                for (j = 0; j < this.m; j++)
                {
                    if (this.matrix[i, j] != 0)
                    {
                        i++;
                        j = -1;
                    }

                    if (i == this.n) break;
                }

                if (j == this.m) return true; //this means that there is a null row in matrix
            }

            for (j = 0; j < this.m; j++)
            {
                for (i = 0; i < this.n; i++)
                {
                    if (this.matrix[i, j] != 0)
                    {
                        j++;
                        i = -1;
                    }

                    if (j == this.m) break;
                }

                if (i == this.n) return true; // means that there is a null column in matrix
            }
            return false;
        }

        /// <summary>
        /// multiplies the matrix with the given value
        /// </summary>
        /// <param name="x">the value to multiply matrix by</param>
        public Matrix Scalar(double value)
        {
            Matrix prod = new Matrix(this.n, this.m);
            for (int i = 0; i < this.n; i++)
            {
                for (int j = 0; j < this.m; j++)
                {
                    prod.matrix[i,j] = this.matrix[i, j] * value;
                }
            }
            return prod;
        }



        //calculates the Inverse of teh matrix if it exists
        /// <summary>
        /// Calculates the inverse of the given matrix
        /// </summary>
        /// <returns>returns the inverse matrix if it exists</returns>
        public Matrix Inverse()
        {
            if (this.n != this.m)
            {
                Console.WriteLine("The inverse for this marix does'nt exit");
                return null;
            }

            if (this.Determinant() == 0)
            {
                Console.WriteLine("Determinant of the matrix is 0.  The inverse for this marix does'nt exit.");
                return null;
            }

            Matrix identity = new Matrix(this.n, 2 * this.n);
            Matrix b = new Matrix(this.n, this.n);
            int i, j, k;
            for (i = 0; i < this.n; i++)
            {
                for (j = 0; j < this.n; j++)
                {
                    identity.matrix[i, j] = this.matrix[i, j];
                }
                for (k = this.n; k < 2 * this.n; k++)
                {
                    if (i == k - this.m)
                        identity.matrix[i, k] = 1;
                    else identity.matrix[i, k] = 0;
                }

            }
            double elemToDivide;
            for (i = 0; i < this.n; i++)
            {
                for (j = i + 1; j < this.n; j++)
                {
                    if (identity.matrix[i, i] == 0)
                    {
                        for (int t = i + 1; t < this.n; t++)
                        {
                            if (identity.matrix[t, i] != 0)
                            {
                                double temp;
                                for (int l = 0; l < this.m; l++)
                                {
                                    temp = this.matrix[i, l];
                                    this.matrix[i, l] = this.matrix[t, l];
                                    this.matrix[l, t] = temp;
                                }
                                break;
                            }
                        }
                    }
                    elemToDivide = identity.matrix[j, i] / identity.matrix[i, i];
                    for (k = i; k < identity.m; k++)
                        identity.matrix[j, k] = identity.matrix[j, k] - elemToDivide * identity.matrix[i, k];
                }
            }

            for (i = this.n - 1; i > 0; i--)
            {
                for (j = i - 1; j >= 0; j--)
                {
                    elemToDivide = identity.matrix[j, i] / identity.matrix[i, i];

                    for (k = 0; k < identity.m; k++)
                        identity.matrix[j, k] = identity.matrix[j, k] - elemToDivide * identity.matrix[i, k];
                }
            }

            for (i = 0; i < identity.n; i++)
            {
                elemToDivide = identity.matrix[i, i];
                for (j = 0; j < identity.m; j++)
                {
                    identity.matrix[i, j] /= elemToDivide;
                }
            }

            Matrix inverseMatrix = new Matrix(this.n, this.m);
            for (i = 0; i < this.n; i++)
            {
                for (j = this.m; j < this.m * 2; j++)
                {
                    inverseMatrix.matrix[i, j - this.m] = identity.matrix[i, j];
                }
            }
            return inverseMatrix;
        }

        /// <summary>
        /// Calculates the transposition of the matrix.
        /// </summary>
        /// <returns>Returns the transpose of a matrix.</returns>
        public Matrix Transpose()
        {
            Matrix t = new Matrix(this.m, this.n);
            for (int i = 0; i < t.n; i++)
            {
                for (int j = 0; j < t.m; j++)
                {
                    t.matrix[i, j] = this.matrix[j, i];
                }
            }
            return t;
        }

        /// <summary>
        /// Check's wheter the matrix is orthogonal or no
        /// </summary>
        /// <returns>returns true if the matrix is orthogonal false otherwise</returns>
        public bool IsOrtogonal()
        {
            if (this.n != this.m)
            {
                Console.WriteLine("This matrix isn't orthogonal");
                return false;
            }
            Matrix a = this * this.Transpose();
            for (int i = 0; i < a.n; i++)
            {
                for (int j = 0; j < a.m; j++)
                {
                    if (j == i && a.matrix[i, j] != 1)
                    {
                        Console.WriteLine("This matrix isn't orthogonal");
                        return false;
                    }
                    else if (a.matrix[i, j] != 0)
                    {

                        Console.WriteLine("This matrix isn't orthogonal");
                        return false;
                    }
                }
            }
            Console.WriteLine("This matrix is orthogonal");
            return true;
        }

        /// <summary>
        /// Translates by the given vector.
        /// </summary>
        /// <param name="vector">the given vector to translate by</param>
        /// <returns>returns the translated vektor</returns>
        public Matrix Translate(Matrix vector)
        {
            if (this.n != vector.n)
            {
                Console.WriteLine("The dimensions of points must be the same.");
                return null;
            }
            Matrix translationMatrix = new Matrix(3, 3);
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    if (i == j)
                        translationMatrix.matrix[i, j] = 1;
                    else translationMatrix.matrix[i, j] = 0;
                }
            }
            for (int i = 0; i < vector.n; i++)
            {
                translationMatrix.matrix[i, vector.n - 1] = vector.matrix[i, 0];
            }

            translationMatrix.matrix[vector.n, vector.n] = 1;

            Matrix duplicate = new Matrix(this.n + 1, this.m);
            for (int i = 0; i < this.n; i++)
                duplicate.matrix[i, 0] = this.matrix[i, 0];
            duplicate.matrix[this.n, 0] = 1;

            Matrix translatedVectorEX = translationMatrix * duplicate;
            Matrix translatedVector = new Matrix(translatedVectorEX.n - 1, translatedVectorEX.m);
            for (int i = 0; i < translatedVectorEX.n - 1; i++)
                translatedVector.matrix[i, 0] = translatedVectorEX.matrix[i, 0];

            return translatedVector;

        }

        /// <summary>
        /// Rotates(2D) by the given angle 
        /// </summary>
        /// <param name="alfa">the given angle to rotate by</param>
        /// <returns>returns the rotated matrix</returns>
        public Matrix Rotate(double alfa)
        {
            double angle = Math.PI * alfa / 180.0;
            Matrix rot = new Matrix(2, 2);
            rot.matrix[0, 0] = Math.Cos(angle);
            rot.matrix[0, 1] = Math.Sin(angle);
            rot.matrix[1, 0] = -Math.Sin(angle);
            rot.matrix[1, 1] = Math.Cos(angle);
            rot.Print();
            Console.WriteLine("its here");
            (rot * this).Print();
            return (this * rot);
        }

        /// <summary>
        /// Scaling(2D) in x, y direction by a factor fact.
        /// </summary>
        /// <returns>returns scale</returns>
        public Matrix Scale()
        {
            Matrix scale = new Matrix(this.n, this.m);
            double val;
            for (int i = 0; i < this.m - 1; i++)
            {
                val = 0;
                for (int j = 0; j < this.n; j++)
                {
                    val += this.matrix[i, j] * this.matrix[i, j];
                }

                scale.matrix[i, i] = Math.Sqrt(val);
            }

            for (int i = 0; i < this.n; i++)
            {
                for (int j = 0; j < this.m; j++)
                {
                    if (i != j) scale.matrix[i, j] = 0;
                }
            }

            scale.matrix[this.n - 1, this.m - 1] = 1;
            return scale;
        }

        /// <summary>
        /// find's the smallest element in matrix
        /// </summary>
        /// <returns>returns the smallest value in the matrix</returns>
        public double Smallest()
        {
            double sm = this.matrix[0, 0];
            for (int i = 0; i < this.n; i++)
            {
                for (int j = 0; j < this.m; j++)
                {
                    if (this.matrix[i, j] < sm)
                        sm = this.matrix[i, j];
                }
            }
            return sm;
        }

        /// <summary>
        /// find's the largest value in the matrix
        /// </summary>
        /// <returns>returns the largest value</returns>
        public double Largest()
        {
            double lg = this.matrix[0, 0];
            for (int i = 0; i < this.n; i++)
            {
                for (int j = 0; j < this.m; j++)
                {
                    if (this.matrix[i, j] > lg)
                        lg = this.matrix[i, j];
                }
            }
            return lg;
        }

        /// <summary>
        /// print's the nxm matrix in the format bellow
        /// m[0,1] m[0,2] ...... m[0,m-1]
        /// m[1,1] m[1,2] ...... m[1,m-1]
        /// ...........................
        /// ...........................
        /// m[n-1,1] m[n-1,2] ...... m[n-1,m-1]
        /// </summary>
        public void Print()
        {
            Console.WriteLine("---------------------");
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    Console.Write(" {0}", this.matrix[i, j]);
                }
                Console.WriteLine("\n");
            }
            Console.WriteLine("---------------------");
        }
    }
}

