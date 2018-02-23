using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
namespace MatrixOp
{
    class Program
    {
        static void Main(string[] args)
        {
            var a = new Matrix(3,3);
            var b = new Matrix(4, 3);
            a.Print();
            b.Print();
            (a+b).Print();
            (a*b).Print();
            (a.Scalar(5)).Print();
            (a.Inverse()).Print();
            a.IsOrtogonal();
            Console.WriteLine("the smallest is {0}",a.Smallest());
            Console.WriteLine("the largest is {0}",a.Largest());
        }
    }
}
