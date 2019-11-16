#pragma once

#include <stdint.h>
#include <Eigen.h>
#include <Eigen/LU>
#include <unsupported/Eigen/MatrixFunctions>

class Magnetometer_Calibration
{
public:
    template <typename Derived>
    void calibrate(const Eigen::ArrayBase<Derived> &mx,
                   const Eigen::ArrayBase<Derived> &my,
                   const Eigen::ArrayBase<Derived> &mz,
                   const Eigen::ArrayBase<Derived> &roll,
                   const Eigen::ArrayBase<Derived> &pitch)
    {
        EIGEN_STATIC_ASSERT_FIXED_SIZE(Derived);
        EIGEN_STATIC_ASSERT(Derived::ColsAtCompileTime == 1,
                            Eigen::THIS_METHOD_IS_ONLY_FOR_MATRICES_OF_A_SPECIFIC_SIZE);

        // 1. fit ellipsoid

        // from https: //de.mathworks.com/matlabcentral/fileexchange/24693-ellipsoid-fit
        /*
        Copyright (c) 2015, Yury Petrov
        All rights reserved.

        Redistribution and use in source and binary forms, with or without
        modification, are permitted provided that the following conditions are
        met:

            * Redistributions of source code must retain the above copyright
            notice, this list of conditions and the following disclaimer.
            * Redistributions in binary form must reproduce the above copyright
            notice, this list of conditions and the following disclaimer in
            the documentation and/or other materials provided with the distribution

        THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
        AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
        IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
        ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
        LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
        CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
        SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
        INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
        CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
        ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
        POSSIBILITY OF SUCH DAMAGE.
        */

        Eigen::Matrix<float, 9, Eigen::Dynamic> D;
        D.resize(9, mx.rows());
        D.row(0) = (mx * mx + my * my - 2 * mz * mz).matrix().transpose();
        D.row(1) = (mx * mx + mz * mz - 2 * my * my).matrix().transpose();
        D.row(2) = (2 * mx * my).matrix().transpose();
        D.row(3) = (2 * mx * mz).matrix().transpose();
        D.row(4) = (2 * mz * my).matrix().transpose();
        D.row(5) = (2 * mx).matrix().transpose();
        D.row(6) = (2 * my).matrix().transpose();
        D.row(7) = (2 * mz).matrix().transpose();
        D.row(8).setOnes();

        Eigen::Matrix<float, Eigen::Dynamic, 1> d2;
        d2.resize(mx.rows(), 1);
        d2 = (mx * mx + my * my + mz * mz).matrix();

        Eigen::Matrix<float, 9, 1> u;
        u = (D * D.transpose()).inverse() * (D * d2);

        Eigen::Matrix<float, 10, 1> v;
        v(0) = u(0) + u(1) - 1;
        v(1) = u(0) - 2 * u(1) - 1;
        v(2) = u(1) - 2 * u(0) - 1;
        v.segment(3, 7) = u.segment(2, 7);

        Eigen::Matrix<float, 4, 4> A;
        A << v(0), v(3), v(4), v(6),
            v(3), v(1), v(5), v(7),
            v(4), v(5), v(2), v(8),
            v(6), v(7), v(8), v(9);

        // V (hard iron) = center of ellipsoid
        V = (-A.topLeftCorner<3, 3>()).inverse() * v.segment<3>(6);

        Eigen::Matrix<float, 4, 4> T;
        T.setIdentity();
        T.row(3).segment(0, 3) = V.transpose();

        Eigen::Matrix<float, 4, 4> R;
        R = T * A * T.transpose();

        Eigen::Matrix<float, 3, 3> M;
        M = R.topLeftCorner<3, 3>().array() / (-R(3, 3));

        float scale = R(3, 3);

        // 2. compute Winv (soft iron)
        Winv = M.sqrt();

        // 3. compute B (geomagnetic field strength)
        B = sqrt(scale);

        // 4. compute inclination
        Eigen::Matrix<float, 3, 3> RxT;
        RxT.setIdentity();
        Eigen::Matrix<float, 3, 3> RyT;
        RyT.setIdentity();
        Eigen::Matrix<float, 3, 1> m;
        Eigen::Matrix<float, 3, 1> Bf;
        inclination = 0;
        for (uint16_t i = 0; i < mx.rows(); i++)
        {
            float croll = cos(roll(i));
            float sroll = sin(roll(i));
            float cpitch = cos(pitch(i));
            float spitch = sin(pitch(i));

            RxT(1, 1) = croll;
            RxT(1, 2) = sroll;
            RxT(2, 1) = -sroll;
            RxT(2, 2) = croll;

            RyT(0, 0) = cpitch;
            RyT(0, 2) = -spitch;
            RyT(2, 0) = spitch;
            RyT(2, 2) = cpitch;

            m(0) = mx(i);
            m(1) = my(i);
            m(2) = mz(i);

            Bf = RxT * RyT * correct(m);
            float yaw = -atan2f(-Bf(1), Bf(0));

            inclination += atanf((-Bf(2) * cos(yaw) / Bf(0))) / mx.rows();
        }
    };

    Eigen::Matrix<float, 3, 1> correct(Eigen::Matrix<float, 3, 1> &m) const
    {
        return Winv * (m - V);
    };

    Eigen::Matrix<float, 3, 3> getWinv() const
    {
        return Winv;
    };

    Eigen::Matrix<float, 3, 3> getW() const
    {
        return Winv.inverse();
    };

    Eigen::Matrix<float, 3, 1> getV() const
    {
        return V;
    };

    float getB() const
    {
        return B;
    };

    float getInclination() const
    {
        return inclination;
    };

protected:
    Eigen::Matrix<float, 3, 3> Winv;
    Eigen::Matrix<float, 3, 1> V;
    float B;
    float inclination;
};