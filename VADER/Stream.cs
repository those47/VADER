using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Intrinsics.Arm;
using System.Text;
using System.Threading.Tasks;

namespace VADER
{
    internal class Stream
    {
        private ulong n;
        private double dt;
        private double dx;

        private double A;
        private double AS;
        private double D;
        private double Q;
        private double qLIN;
        private double alpha;
        private double ChatS;
        private double Kd;
        private double lambda;
        private double lambdaS;
        private double lambdahat;
        private double lambdahatS;
        private double rho;
        private double CL;

        private double[,] C;
        private double[,] CS;
        private double[,] Csed;

        private double gamma;
        private double G;
        private double E;
        private double F;

        public Stream(ulong streamLength, double totalDuration, double dt, double dx, double[] upstreamBoundaryData, double[] parameters)
        {
            ulong numberOfSegments = (ulong)(streamLength / dx);

            this.n = numberOfSegments;
            this.dt = dt;
            this.dx = dx;

            this.A = parameters[0];
            this.AS = parameters[1];
            this.CL = parameters[2];
            this.D = parameters[3];
            this.Q = parameters[4];
            this.qLIN = parameters[5];
            this.alpha = parameters[6];
            this.ChatS = parameters[7];
            this.Kd = parameters[8];
            this.lambda = parameters[9];
            this.lambdaS = parameters[10];
            this.lambdahat = parameters[11];
            this.lambdahatS = parameters[12];
            this.rho = parameters[13];

            this.gamma = alpha * dt * A / AS;
            this.G = dt / (2 * A * dx) * (Q / 2 - A * D / dx);
            this.E = -dt / (2 * A * dx) * (Q / 2 + A * D / dx);
            this.F = 1 + dt / 2 * ((A * D + A * D) / (A * dx * dx) + qLIN / A
                    + alpha * (1 - alpha * dt * A / AS / (2 + alpha * dt * A / AS + dt * lambdahatS + dt * lambdaS))
                    + rho * lambdahat * Kd * (1 - dt * lambdahat / (2 + dt * lambdahat)) + lambda);

            this.C = new double[numberOfSegments, (ulong)(totalDuration / dt)];
            this.CS = new double[numberOfSegments, (ulong)(totalDuration / dt)];
            this.Csed = new double[numberOfSegments, (ulong)(totalDuration / dt)];
            for (int i = 1; i < C.GetLength(0); i++)
            {
                C[i, 0] = 0;
                CS[i, 0] = 0;
                Csed[i, 0] = 0;
            }

            for (int j = 0; j < C.GetLength(1) && j < upstreamBoundaryData.Length; j++)
            {
                C[0, j] = upstreamBoundaryData[j];
                CS[0, j] = 0;
                Csed[0, j] = 0;
            }
            for (int j = upstreamBoundaryData.Length; j < C.GetLength(1); j++)
            {
                C[0, j] = 0;
                CS[0, j] = 0;
                Csed[0, j] = 0;
            }
        }

        private double R(int i, int j)
        {
            return C[i, j] + dt / 2 * (G + qLIN / A * CL
                    + alpha * ((2 - gamma - dt * lambdahatS - dt * lambdaS)* CS[i, j] + gamma * C[i, j] + 2 * dt * lambdahatS * ChatS)
                        / (2 + gamma + dt * lambdahatS + dt * lambdaS))
                    + rho * lambdahat * ((2 - dt * lambdahat) * Csed[i, j] + dt * lambdahat * Kd * C[i, j]) / (2 + dt * lambdahat);
        }

        public List<double[,]> SimulateAdvectionAndDispersion()
        {
            PropagateTime();
            return new List<double[,]> { C, CS, Csed };
        }

        private void PropagateTime()
        {
            for (int j = 0; j < C.GetLength(1) - 1; j++)
            {
                IncrementTimeForward(j);
            }
        }
        // Propagates C, CS, and Csed from time j to time j+1 -- See https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
        
        private void IncrementTimeForward(int j)
        {
            // Set up the RHS of the Thomas algorithm matrix equation
            double[] d = new double[C.GetLength(0) - 1];

            // Set up the LHS of the Thomas algorithm matrix equation
            double[] a = new double[C.GetLength(0) - 1];
            double[] b = new double[C.GetLength(0) - 1];
            double[] c = new double[C.GetLength(0) - 1];


            // Populate the values of the matrix element and vector arrays
            for (int i = 0; i < C.GetLength(0) - 1; i++) 
            {
                d[i] = R(i, j);
                a[i] = E;
                b[i] = F;
                c[i] = G;
            }
            // This line enforces the upstream boundary condition of the given upstream concentration
            d[0] = R(1, j) - E * C[0, j + 1];
            //This line enforces the downstream boundary condition of no change in concentration along the stream
            b[C.GetLength(0) - 2] = F + G;

            // Perform the Thomas' agorithm
            
            // Modify the coefficients
            double[] cp = new double[C.GetLength(0) - 2];
            cp[0] = c[0] / b[0];
            double[] dp = new double[C.GetLength(0) - 1];
            dp[0] = d[0] / b[0];

            for (int i = 1; i < C.GetLength(0) - 2; i++)
            {
                cp[i] = c[i] / (b[i] - a[i] * cp[i - 1]);
                dp[i] = (d[i] - a[i] * dp[i - 1]) / (b[i] - a[i] * cp[i - 1]);
            }
            dp[C.GetLength(0) - 2] = (d[C.GetLength(0) - 2] - a[C.GetLength(0) - 2] * dp[C.GetLength(0) - 2 - 1]) / (b[C.GetLength(0) - 2] - a[C.GetLength(0) - 2] * cp[C.GetLength(0) - 2 - 1]);


            // Perform the backsubstitution to obtain the new concentrations
            C[C.GetLength(0) - 1, j + 1] = dp[C.GetLength(0) - 2];
            CS[C.GetLength(0) - 1, j + 1] = ((2 - gamma - dt * lambdahatS - dt * lambdaS) * CS[C.GetLength(0) - 1, j] + gamma * C[C.GetLength(0) - 1, j] + gamma * C[C.GetLength(0) - 1, j + 1] + 2 * dt * lambdahatS * ChatS)
                                            / (2 + gamma + dt * lambdahatS + dt * lambdaS);
            Csed[C.GetLength(0) - 1, j + 1] = ((2 - dt * lambdahat) * Csed[C.GetLength(0) - 1, j] + dt * lambdahat * Kd * (C[C.GetLength(0) - 1, j] + C[C.GetLength(0) - 1, j+1])) / (2 + dt * lambdahat);
            for (int i = C.GetLength(0) - 2; i >= 1; i--)
            {
                C[i, j + 1] = dp[i - 1] - cp[i - 1] * C[i + 1, j + 1];
                CS[i, j + 1] = ((2 - gamma - dt * lambdahatS - dt * lambdaS) * CS[i, j] + gamma * C[i, j] + gamma * C[i, j + 1] + 2 * dt * lambdahatS * ChatS)
                                                / (2 + gamma + dt * lambdahatS + dt * lambdaS);
                Csed[i, j + 1] = ((2 - dt * lambdahat) * Csed[i, j] + dt * lambdahat * Kd * (C[i, j] + C[i, j + 1])) / (2 + dt * lambdahat);
            }
        }


    }
}
