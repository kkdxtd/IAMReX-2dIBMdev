// SPDX-FileCopyrightText: 2023 - 2025 Yadong Zeng<zdsjtu@gmail.com> & ZhuXu Li<1246206018@qq.com>
//
// SPDX-License-Identifier: BSD-3-Clause

#include <AMReX_Math.H>
#include <AMReX_Print.H>
#include <DiffusedIB.H>

#include <AMReX_ParmParse.H>
#include <AMReX_TagBox.H>
#include <AMReX_Utility.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_FillPatchUtil.H>
#include <iamr_constants.H>

#include <filesystem>
#include <sstream>
namespace fs = std::filesystem;

#define GHOST_CELLS 2

using namespace amrex;

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                     global variable                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#define LOCAL_LEVEL 0
#if (AMREX_SPACEDIM == 3)
const Vector<std::string> direction_str{"X","Y","Z"};
namespace ParticleProperties{
    Vector<Real> _x{}, _y{}, _z{}, _rho{};
    Vector<Real> Vx{}, Vy{}, Vz{};
    Vector<Real> Ox{}, Oy{}, Oz{};
    Vector<Real> _radius;
    Real rd{0.0};
    Vector<int> TLX{}, TLY{},TLZ{},RLX{},RLY{},RLZ{};
    int euler_finest_level{0};
    int euler_velocity_index{0};
    int euler_force_index{0};
    Real euler_fluid_rho{0.0};
    int verbose{0};
    int loop_ns{2};
    int loop_solid{1};
    int Uhlmann{0};

    Vector<Real> GLO, GHI;
    int start_step{-1};
    int collision_model{0};

    int write_freq{1};
    bool init_particle_from_file{false};

    GpuArray<Real, 3> plo{0.0,0.0,0.0}, phi{0.0,0.0,0.0}, dx{0.0, 0.0, 0.0};
}
#endif
#if (AMREX_SPACEDIM == 2)
namespace FiberProperties{
    Vector<Real> x0{}, y0{};                    // 细丝起点坐标
    Vector<Real> amplitude{};                   // 振幅A
    Vector<Real> period{};                      // 周期T
    Vector<Real> length{};                      // 长度L
    Vector<Real> wave_number{};                 // 波数n1
    Vector<Real> phase{};                       // 初始相位phi1
    Vector<int> num_marker{};                   // 拉格朗日点数
    int euler_finest_level{0};                  // 最细网格等级
    int euler_velocity_index{0};                // 欧拉速度场索引
    int euler_force_index{0};                   // 欧拉力场索引
    Real euler_fluid_rho{0.0};                  // 流体密度
    int verbose{0};                             // 调试输出
    int loop_ns{2};                             // 欧拉-拉格朗日交互迭代次数
    int start_step{-1};                         // 开始步数
    int write_freq{1};                          // 输出频率
    Vector<Real> GLO, GHI;                      // 几何边界
    GpuArray<Real, 2> plo{0.0,0.0}, phi{0.0,0.0}, dx{0.0, 0.0};
}
#endif

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                     other function                            */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#if (AMREX_SPACEDIM == 3)
void nodal_phi_to_pvf(MultiFab& pvf, const MultiFab& phi_nodal)
{

    // amrex::Print() << "In the nodal_phi_to_pvf\n";

    pvf.setVal(0.0);

    // Only set the valid cells of pvf
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(pvf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto const& pvffab   = pvf.array(mfi);
        auto const& pnfab = phi_nodal.array(mfi);
        amrex::ParallelFor(bx, [pvffab, pnfab]
        AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            Real num = 0.0;
            for(int kk=k; kk<=k+1; kk++) {
                for(int jj=j; jj<=j+1; jj++) {
                    for(int ii=i; ii<=i+1; ii++) {
                        num += (-pnfab(ii,jj,kk)) * nodal_phi_to_heavi(-pnfab(ii,jj,kk));
                    }
                }
            }
            Real deo = 0.0;
            for(int kk=k; kk<=k+1; kk++) {
                for(int jj=j; jj<=j+1; jj++) {
                    for(int ii=i; ii<=i+1; ii++) {
                        deo += std::abs(pnfab(ii,jj,kk));
                    }
                }
            }
            pvffab(i,j,k) = num / (deo + 1.e-12);
        });
    }

}

void calculate_phi_nodal(MultiFab& phi_nodal, kernel& current_kernel)
{
    phi_nodal.setVal(0.0);

    amrex::Real Xp = current_kernel.location[0];
    amrex::Real Yp = current_kernel.location[1];
    amrex::Real Zp = current_kernel.location[2];
    amrex::Real Rp = current_kernel.radius;

    // Only set the valid cells of phi_nodal
    for (MFIter mfi(phi_nodal,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto const& pnfab = phi_nodal.array(mfi);
        auto dx = ParticleProperties::dx;
        auto plo = ParticleProperties::plo;
        amrex::ParallelFor(bx, [=]
            AMREX_GPU_DEVICE(int i, int j, int k) noexcept
            {
                Real Xn = i * dx[0] + plo[0];
                Real Yn = j * dx[1] + plo[1];
                Real Zn = k * dx[2] + plo[2];

                pnfab(i,j,k) = std::sqrt( (Xn - Xp)*(Xn - Xp)
                        + (Yn - Yp)*(Yn - Yp)  + (Zn - Zp)*(Zn - Zp)) - Rp;
                pnfab(i,j,k) = pnfab(i,j,k) / Rp;

            }
        );
    }
}

// May use ParReduce later, https://amrex-codes.github.io/amrex/docs_html/GPU.html#multifab-reductions
void CalculateSumU_cir (RealVect& sum,
                        const MultiFab& E,
                        const MultiFab& pvf,
                        int EulerVelIndex)
{
    auto const& E_data = E.const_arrays();
    auto const& pvf_data = pvf.const_arrays();
    const Real d = Math::powi<3>(ParticleProperties::dx[0]);
    amrex::GpuTuple<Real, Real, Real> tmpSum = ParReduce(TypeList<ReduceOpSum,ReduceOpSum,ReduceOpSum>{}, TypeList<Real, Real, Real>{},E, IntVect{0},
    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> amrex::GpuTuple<Real, Real, Real>{
        auto E_ = E_data[box_no];
        auto pvf_ = pvf_data[box_no];
        return {
            E_(i, j, k, EulerVelIndex    ) * d * pvf_(i,j,k),
            E_(i, j, k, EulerVelIndex + 1) * d * pvf_(i,j,k),
            E_(i, j, k, EulerVelIndex + 2) * d * pvf_(i,j,k)
        };
    });
    sum[0] = amrex::get<0>(tmpSum);
    sum[1] = amrex::get<1>(tmpSum);
    sum[2] = amrex::get<2>(tmpSum);
}

void CalculateSumT_cir (RealVect& sum,
                        const MultiFab& E,
                        const MultiFab& pvf,
                        const RealVect pLoc,
                        int EulerVelIndex)
{
    auto plo = ParticleProperties::plo;
    auto dx = ParticleProperties::dx;

    auto const& E_data = E.const_arrays();
    auto const& pvf_data = pvf.const_arrays();
    const Real d = Math::powi<3>(ParticleProperties::dx[0]);
    amrex::GpuTuple<Real, Real, Real> tmpSum = ParReduce(TypeList<ReduceOpSum,ReduceOpSum,ReduceOpSum>{}, TypeList<Real, Real, Real>{},E, IntVect{0},
    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> amrex::GpuTuple<Real, Real, Real>{
        auto E_ = E_data[box_no];
        auto pvf_ = pvf_data[box_no];

        Real x = plo[0] + i*dx[0] + 0.5*dx[0];
        Real y = plo[1] + j*dx[1] + 0.5*dx[1];
        Real z = plo[2] + k*dx[2] + 0.5*dx[2];

        Real vx = E_(i, j, k, EulerVelIndex    );
        Real vy = E_(i, j, k, EulerVelIndex + 1);
        Real vz = E_(i, j, k, EulerVelIndex + 2);

        RealVect tmp = RealVect(x - pLoc[0], y - pLoc[1], z - pLoc[2]).crossProduct(RealVect(vx, vy, vz));

        return {
            tmp[0] * d * pvf_(i, j, k),
            tmp[1] * d * pvf_(i, j, k),
            tmp[2] * d * pvf_(i, j, k)
        };
    });
    sum[0] = amrex::get<0>(tmpSum);
    sum[1] = amrex::get<1>(tmpSum);
    sum[2] = amrex::get<2>(tmpSum);
}

[[nodiscard]] AMREX_FORCE_INLINE
Real cal_momentum(Real rho, Real radius)
{
    return 8.0 * Math::pi<Real>() * rho * Math::powi<5>(radius) / 15.0;
}
#endif
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void deltaFunction(Real xf, Real xp, Real h, Real& value, DELTA_FUNCTION_TYPE type)
{
    Real rr = amrex::Math::abs(( xf - xp ) / h);

    switch (type) {
    case DELTA_FUNCTION_TYPE::FOUR_POINT_IB:
        if(rr >= 0 && rr < 1 ){
            value = 1.0 / 8.0 * ( 3.0 - 2.0 * rr + std::sqrt( 1.0 + 4.0 * rr - 4 * Math::powi<2>(rr))) / h;
        }else if (rr >= 1 && rr < 2) {
            value = 1.0 / 8.0 * ( 5.0 - 2.0 * rr - std::sqrt( -7.0 + 12.0 * rr - 4 * Math::powi<2>(rr))) / h;
        }else {
            value = 0;
        }
        break;
    case DELTA_FUNCTION_TYPE::THREE_POINT_IB:
        if(rr >= 0.5 && rr < 1.5){
            value = 1.0 / 6.0 * ( 5.0 - 3.0 * rr - std::sqrt( - 3.0 * Math::powi<2>( 1 - rr) + 1.0 )) / h;
        }else if (rr >= 0 && rr < 0.5) {
            value = 1.0 / 3.0 * ( 1.0 + std::sqrt( 1.0 - 3 * Math::powi<2>(rr))) / h;
        }else {
            value = 0;
        }
        break;
    default:
        break;
    }
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                    mParticle/mFiber member function                  */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
//loop all particels
#if (AMREX_SPACEDIM == 3)
void mParticle::InteractWithEuler(MultiFab &EulerVel, 
                                  MultiFab &EulerForce, 
                                  Real dt,
                                  DELTA_FUNCTION_TYPE type)
{
    if (verbose) amrex::Print() << "[Particle] mParticle::InteractWithEuler\n";
    // clear time , start record
    spend_time = 0;
    auto InteractWithEulerStart = ParallelDescriptor::second();

    MultiFab EulerForceTmp(EulerForce.boxArray(), EulerForce.DistributionMap(), 3, EulerForce.nGrow());
    //clean preStep's IB_porperties 
    for(auto& kernel : particle_kernels) {
        kernel.ib_force.scale(0.0);
        kernel.ib_moment.scale(0.0);
    }

    //for 1 -> Ns
    int loop = ParticleProperties::loop_ns;
    BL_ASSERT(loop > 0);
    while(loop > 0){
        if(verbose) amrex::Print() << "[Particle] Ns loop index : " << loop << "\n";
        
        EulerForce.setVal(0.0);

        for(kernel& kernel : particle_kernels){
            InitialWithLargrangianPoints(kernel); // Initialize markers for a specific particle
            ResetLargrangianPoints();
            EulerForceTmp.setVal(0.0);
            auto ib_force = kernel.ib_force;
            auto ib_moment = kernel.ib_moment;
            kernel.ib_force.scale(0.0); // clear kernel ib_force
            kernel.ib_moment.scale(0.0); // clear kernel ib_moment

            VelocityInterpolation(EulerVel, type);
            ComputeLagrangianForce(dt, kernel);
            ForceSpreading(EulerForceTmp, kernel, type);
            MultiFab::Add(EulerForce, EulerForceTmp, 0, 0, 3, EulerForce.nGrow());

            kernel.ib_force += ib_force;
            kernel.ib_moment += ib_moment;
        }
        VelocityCorrection(EulerVel, EulerForce, dt);
        loop--;
    }
    spend_time += ParallelDescriptor::second() - InteractWithEulerStart;
}
#endif
#if (AMREX_SPACEDIM == 2)
void mFiber::InteractWithEuler(MultiFab &EulerVel, MultiFab &EulerForce, Real dt, DELTA_FUNCTION_TYPE type)
{
    if (verbose) amrex::Print() << "[Fiber] mFiber::InteractWithEuler\n";
    
    // 清零时间记录
    spend_time = 0;
    auto InteractWithEulerStart = ParallelDescriptor::second();

    MultiFab EulerForceTmp(EulerForce.boxArray(), EulerForce.DistributionMap(), 2, EulerForce.nGrow());

    
    // 清零所有细丝的IB力
    for(auto& fiber_kernel : fiber_kernels) {
        fiber_kernel.ib_force.scale(0.0);
    }

    // 欧拉-拉格朗日交互主循环
    int loop = FiberProperties::loop_ns;
    BL_ASSERT(loop > 0);
    while(loop > 0){
        if(verbose) amrex::Print() << "[Fiber] Ns loop index : " << loop << "\n";
        
        EulerForce.setVal(0.0);

        for(fiber_kernel& fiber_kernel : fiber_kernels){
            // 初始化拉格朗日点位置
            InitialWithLargrangianPoints(fiber_kernel);
            ResetLargrangianPoints();
            EulerForceTmp.setVal(0.0);
            
            // 保存旧的IB力
            auto ib_force = fiber_kernel.ib_force;
            fiber_kernel.ib_force.scale(0.0);

            // IBM核心步骤
            VelocityInterpolation(EulerVel, type);
            ComputeLagrangianForce(dt, fiber_kernel);
            ForceSpreading(EulerForceTmp, fiber_kernel, type);
            
            // 累加所有细丝的力到总欧拉力场
            MultiFab::Add(EulerForce, EulerForceTmp, 0, 0, 2, EulerForce.nGrow());

            // 恢复IB力
            fiber_kernel.ib_force += ib_force;
        }
        
        // 速度修正
        VelocityCorrection(EulerVel, EulerForce, dt);
        loop--;
    }
    
    spend_time += ParallelDescriptor::second() - InteractWithEulerStart;
}
#endif
#if (AMREX_SPACEDIM == 3)
void mParticle::InitParticles(const Vector<Real>& x,
                              const Vector<Real>& y,
                              const Vector<Real>& z,
                              const Vector<Real>& rho_s,
                              const Vector<Real>& Vx,
                              const Vector<Real>& Vy,
                              const Vector<Real>& Vz,
                              const Vector<Real>& Ox,
                              const Vector<Real>& Oy,
                              const Vector<Real>& Oz,
                              const Vector<int>& TLXt,
                              const Vector<int>& TLYt,
                              const Vector<int>& TLZt,
                              const Vector<int>& RLXt,
                              const Vector<int>& RLYt,
                              const Vector<int>& RLZt,
                              const Vector<Real>& radius,
                              Real h,
                              Real gravity,
                              int _verbose)
{
    verbose = _verbose;
    if (verbose) amrex::Print() << "[Particle] mParticle::InitParticles\n";

    m_gravity[2] = gravity;

    //pre judge
    if(!((x.size() == y.size()) && (x.size() == z.size()))){
        Print() << "particle's position container are all different size";
        return;
    }

    //all the particles have different radius
    for(int index = 0; index < x.size(); index++){
        int real_index;
        // if initial with input file, initialized by [0] data
        if(ParticleProperties::init_particle_from_file){
            real_index = 0;
        }else{
            real_index = index;
        }

        kernel mKernel;
        mKernel.id = index + 1;
        mKernel.location[0] = x[index];
        mKernel.location[1] = y[index];
        mKernel.location[2] = z[index];
        mKernel.velocity[0] = Vx[real_index];
        mKernel.velocity[1] = Vy[real_index];
        mKernel.velocity[2] = Vz[real_index];
        mKernel.omega[0] = Ox[real_index];
        mKernel.omega[1] = Oy[real_index];
        mKernel.omega[2] = Oz[real_index];

        // use current state to initialize old state
        mKernel.location_old = mKernel.location;
        mKernel.velocity_old = mKernel.velocity;
        mKernel.omega_old = mKernel.omega;

        mKernel.TL[0] = TLXt[real_index];
        mKernel.TL[1] = TLYt[real_index];
        mKernel.TL[2] = TLZt[real_index];
        mKernel.RL[0] = RLXt[real_index];
        mKernel.RL[1] = RLYt[real_index];
        mKernel.RL[2] = RLZt[real_index];
        mKernel.rho = rho_s[real_index];
        mKernel.radius = radius[real_index];
        mKernel.Vp = Math::pi<Real>() * 4 / 3 * Math::powi<3>(radius[real_index]);

        //int Ml = static_cast<int>( Math::pi<Real>() / 3 * (12 * Math::powi<2>(mKernel.radius / h)));
        //Real dv = Math::pi<Real>() * h / 3 / Ml * (12 * mKernel.radius * mKernel.radius + h * h);
        int Ml = static_cast<int>((amrex::Math::powi<3>(mKernel.radius - (ParticleProperties::rd - 0.5) * h)
               - amrex::Math::powi<3>(mKernel.radius - (ParticleProperties::rd + 0.5) * h))/(3.*h*h*h/4./Math::pi<Real>()));
        Real dv = (amrex::Math::powi<3>(mKernel.radius - (ParticleProperties::rd - 0.5) * h)
               - amrex::Math::powi<3>(mKernel.radius - (ParticleProperties::rd + 0.5) * h))/(3.*Ml/4./Math::pi<Real>());
        mKernel.ml = Ml;
        mKernel.dv = dv;
        if( Ml > max_largrangian_num ) max_largrangian_num = Ml;

        Real phiK = 0;
        for(int marker_index = 0; marker_index < Ml; marker_index++){
            Real Hk = -1.0 + 2.0 * (marker_index) / ( Ml - 1.0);
            Real thetaK = std::acos(Hk);    
            if(marker_index == 0 || marker_index == (Ml - 1)){
                phiK = 0;
            }else {
                phiK = std::fmod( phiK + 3.809 / std::sqrt(Ml) / std::sqrt( 1 - Math::powi<2>(Hk)) , 2 * Math::pi<Real>());
            }
            mKernel.phiK.push_back(phiK);
            mKernel.thetaK.push_back(thetaK);
        }

        particle_kernels.emplace_back(mKernel);

        if (verbose) amrex::Print() << "h: " << h << ", Ml: " << Ml << ", D: " << Math::powi<3>(h) << " gravity : " << gravity << "\n"
                                    << "Kernel : " << index << ": Location (" << x[index] << ", " << y[index] << ", " << z[index] 
                                    << "), Velocity : (" << mKernel.velocity[0] << ", " << mKernel.velocity[1] << ", "<< mKernel.velocity[2] 
                                    << "), Radius: " << mKernel.radius << ", Ml: " << Ml << ", dv: " << dv << ", Rho: " << mKernel.rho << "\n";
    }
    //collision box generate
    m_Collision.SetGeometry(RealVect(ParticleProperties::GLO), RealVect(ParticleProperties::GHI),particle_kernels[0].radius, h);
}
#endif
#if (AMREX_SPACEDIM == 2)
void mFiber::InitFibers(const Vector<Real>& x0,
                         const Vector<Real>& y0,
                         const Vector<Real>& amplitude,
                         const Vector<Real>& period,
                         const Vector<Real>& length,
                         const Vector<Real>& wave_number,
                         const Vector<Real>& phase,
                         const Vector<int>& num_marker,
                         Real h,
                         int _verbose)
{
    verbose = _verbose;
    if (verbose) amrex::Print() << "[Fiber] mFiber::InitFibers\n";
    
    for(int index = 0; index < x0.size(); index++){
        fiber_kernel mKernel;
        mKernel.id = index + 1;
        mKernel.amplitude = amplitude[index];
        mKernel.period = period[index];
        mKernel.length = length[index];
        mKernel.wave_number = wave_number[index];
        mKernel.phase = phase[index];
        mKernel.num_marker = num_marker[index];

        mKernel.y0_base = y0[index];  //保存基线

        //计算拉格朗日点数，确保两端端点都有点，间距为网格长度h
        int num_points = static_cast<int>(std::ceil(length[index] / h)) + 1;
        mKernel.num_marker = num_points;
        
        // 初始化拉格朗日点位置（使用计算得到的 mKernel.num_marker 保持一致）
        for (int i = 0; i < mKernel.num_marker; i++) {
            Real x_pos = x0[index] + i * h;
            Real y_pos = y0[index] + amplitude[index] *
                         sin(2 * M_PI * (-wave_number[index] * x_pos / length[index]) + phase[index]);

            mKernel.x_marker.push_back(x_pos);
            mKernel.y_marker.push_back(y_pos);
            mKernel.velocity_marker.push_back(RealVect(0.0, 0.0));
        }
        
        // 计算面积 da，仿照三维颗粒中 dv 的计算方式
        Real da = h * h;// 每个拉格朗日点对应的面积，假设为 h^2
        mKernel.da = da;

        //mKernel.x_marker_d.assign(mKernel.x_marker.begin(), mKernel.x_marker.end());
        //mKernel.y_marker_d.assign(mKernel.y_marker.begin(), mKernel.y_marker.end());
        //mKernel.velocity_marker_d.assign(mKernel.velocity_marker.begin(), mKernel.velocity_marker.end());
        // 同步到 device
        mKernel.x_marker_d.resize(mKernel.x_marker.size());
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                              mKernel.x_marker.begin(), mKernel.x_marker.end(),
                              mKernel.x_marker_d.begin());

        mKernel.y_marker_d.resize(mKernel.y_marker.size());
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                              mKernel.y_marker.begin(), mKernel.y_marker.end(),
                              mKernel.y_marker_d.begin());

        mKernel.velocity_marker_d.resize(mKernel.velocity_marker.size());
        amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                              mKernel.velocity_marker.begin(), mKernel.velocity_marker.end(),
                              mKernel.velocity_marker_d.begin());
        amrex::Gpu::streamSynchronize();

        fiber_kernels.emplace_back(mKernel);
        
        if (verbose) amrex::Print() << "Fiber " << index << ": x0=" << x0[index] 
                                    << ", y0=" << y0[index] << ", amplitude=" << amplitude[index]
                                    << ", period=" << period[index] << ", length=" << length[index]
                                    << ", num_marker=" << num_marker[index] << "\n";
    }
}
#endif
#if (AMREX_SPACEDIM == 3)
void mParticle::InitialWithLargrangianPoints(const kernel& current_kernel){

    if (verbose) amrex::Print() << "mParticle::InitialWithLargrangianPoints\n";
    for(mParIter pti(*mContainer, LOCAL_LEVEL); pti.isValid(); ++pti){
        const Long np = pti.numParticles();
        if(np == 0) continue;
        auto *particles = pti.GetArrayOfStructs().data();

        const auto location = current_kernel.location;
        const auto radius = current_kernel.radius;
        const auto* phiK = current_kernel.phiK.dataPtr();
        const auto* thetaK = current_kernel.thetaK.dataPtr();

        amrex::ParallelFor( np, [=]
            AMREX_GPU_DEVICE (int i) noexcept {
                auto id = particles[i].id();
                particles[i].pos(0) = location[0] + radius * std::sin(thetaK[id - 1]) * std::cos(phiK[id - 1]);
                particles[i].pos(1) = location[1] + radius * std::sin(thetaK[id - 1]) * std::sin(phiK[id - 1]);
                particles[i].pos(2) = location[2] + radius * std::cos(thetaK[id - 1]);
            }
        );
    }
    // Redistribute the markers after updating their locations
    mContainer->Redistribute();
    if (verbose) {
        amrex::Print() << "[particle] : particle num :" << mContainer->TotalNumberOfParticles() << "\n";
        mContainer->WriteAsciiFile(amrex::Concatenate("particle", 1));
    }
}
#endif
#if (AMREX_SPACEDIM == 2)
void mFiber::InitialWithLargrangianPoints(const fiber_kernel& kernel)
{
    if (verbose) amrex::Print() << "mFiber::InitialWithLargrangianPoints\n";
    
    for(mFiberIter pti(*mContainer, LOCAL_LEVEL); pti.isValid(); ++pti){
        const Long np = pti.numParticles();
        if(np == 0) continue;
        auto *particles = pti.GetArrayOfStructs().data();

        // 获取细丝的拉格朗日点坐标
        //const auto* x_marker = kernel.x_marker.data();
        //const auto* y_marker = kernel.y_marker.data();
        const auto* x_marker = kernel.x_marker_d.dataPtr();
        const auto* y_marker = kernel.y_marker_d.dataPtr();

        amrex::ParallelFor( np, [=]
            AMREX_GPU_DEVICE (int i) noexcept {
                auto id = particles[i].id();
                // 设置拉格朗日点位置（沿细丝长度方向分布）
                particles[i].pos(0) = x_marker[id - 1];
                particles[i].pos(1) = y_marker[id - 1];
            }
        );
    }
    
    // 重新分布拉格朗日点
    mContainer->Redistribute();
    
    if (verbose) {
        amrex::Print() << "[fiber] : fiber marker num :" << mContainer->TotalNumberOfParticles() << "\n";
        mContainer->WriteAsciiFile(amrex::Concatenate("fiber", 1));
    }

    //if (verbose) {
        // 统计全局总数，而非单个进程本地数
    //    long local_cnt = static_cast<long>(mContainer->NumberOfParticles());
    //    ParallelDescriptor::ReduceLongSum(local_cnt);
    //    amrex::Print() << "[fiber] : fiber marker num (global) :" << local_cnt << "\n";
    //    mContainer->WriteAsciiFile(amrex::Concatenate("fiber", 1));
    //}
}
#endif
#if (AMREX_SPACEDIM == 3)
template <typename P = Particle<numAttri>>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void VelocityInterpolation_cir(int p_iter, P const& p, Real& Up, Real& Vp, Real& Wp,
                               Array4<Real const> const& E, int EulerVIndex,
                               const int *lo, const int *hi, 
                               GpuArray<Real, AMREX_SPACEDIM> const& plo,
                               GpuArray<Real, AMREX_SPACEDIM> const& dx,
                               DELTA_FUNCTION_TYPE type)
{

    //std::cout << "lo " << lo[0] << " " << lo[1] << " "<< lo[2] << "\n";
    //std::cout << "hi " << hi[0] << " " << hi[1] << " "<< hi[2] << "\n";

    const Real d = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);

    const Real lx = (p.pos(0) - plo[0]) / dx[0]; // x
    const Real ly = (p.pos(1) - plo[1]) / dx[1]; // y
    const Real lz = (p.pos(2) - plo[2]) / dx[2]; // z

    int i = static_cast<int>(Math::floor(lx)); // i
    int j = static_cast<int>(Math::floor(ly)); // j
    int k = static_cast<int>(Math::floor(lz)); // k

    //std::cout << "p_iter " << p_iter << " p.pos(0): " << p.pos(0) << " p.pos(1): " << p.pos(1) << " p.pos(2): " << p.pos(2) << "\n";

    // std::cout << "d: " << d << "\n"
    //         << "lx: " << lx << ", ly: " << ly << ", lz: " << lz << "\n"
    //         << "i: " << i << ", j: " << j << ", k: " << k << std::endl;

    Up = 0;
    Vp = 0;
    Wp = 0;
    //Euler to largrangian
    for(int ii = -2; ii < 3; ii++){
        for(int jj = -2; jj < 3; jj++){
            for(int kk = -2; kk < 3; kk ++){
                Real tU, tV, tW;
                const Real xi = plo[0] + (i + ii) * dx[0] + dx[0]/2;
                const Real yj = plo[1] + (j + jj) * dx[1] + dx[1]/2;
                const Real kz = plo[2] + (k + kk) * dx[2] + dx[2]/2;
                deltaFunction( p.pos(0), xi, dx[0], tU, type);
                deltaFunction( p.pos(1), yj, dx[1], tV, type);
                deltaFunction( p.pos(2), kz, dx[2], tW, type);
                const Real delta_value = tU * tV * tW;
                Up += delta_value * E(i + ii, j + jj, k + kk, EulerVIndex    ) * d;
                Vp += delta_value * E(i + ii, j + jj, k + kk, EulerVIndex + 1) * d;
                Wp += delta_value * E(i + ii, j + jj, k + kk, EulerVIndex + 2) * d;
            }
        }
    }
}
#endif
#if (AMREX_SPACEDIM == 2)
template <typename P = Particle<numAttri>>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void VelocityInterpolation_cir(int p_iter, P const& p, 
                               Real& Up, Real& Vp,
                               Array4<Real const> const& E, int EulerVIndex,
                               const int *lo, const int *hi, 
                               GpuArray<Real, AMREX_SPACEDIM> const& plo,
                               GpuArray<Real, AMREX_SPACEDIM> const& dx,
                               DELTA_FUNCTION_TYPE type)
{
    // 体积计算
    const Real d = dx[0] * dx[1];

    // 坐标计算
    const Real lx = (p.pos(0) - plo[0]) / dx[0];
    const Real ly = (p.pos(1) - plo[1]) / dx[1];
    int i = static_cast<int>(Math::floor(lx));
    int j = static_cast<int>(Math::floor(ly));


    // 初始化速度分量
    Up = 0; Vp = 0;


    // 循环
    for(int ii = -2; ii < 3; ii++){
        for(int jj = -2; jj < 3; jj++){
            Real tU, tV;
            const Real xi = plo[0] + (i + ii) * dx[0] + dx[0]/2;
            const Real yj = plo[1] + (j + jj) * dx[1] + dx[1]/2;
            deltaFunction(p.pos(0), xi, dx[0], tU, type);
            deltaFunction(p.pos(1), yj, dx[1], tV, type);
            const Real delta_value = tU * tV;
            Up += delta_value * E(i + ii, j + jj, 0, EulerVIndex    ) * d;
            Vp += delta_value * E(i + ii, j + jj, 0, EulerVIndex + 1) * d;
        }
    }
}
#endif

#if (AMREX_SPACEDIM == 3)
void mParticle::VelocityInterpolation(MultiFab &EulerVel,
                                      DELTA_FUNCTION_TYPE type)//
{
    if (verbose) amrex::Print() << "\tmParticle::VelocityInterpolation\n";

    //amrex::Print() << "euler_finest_level " << euler_finest_level << std::endl;
    const auto& gm = mContainer->GetParGDB()->Geom(LOCAL_LEVEL);
    auto plo = gm.ProbLoArray();
    auto dx = gm.CellSizeArray();
    // attention
    // velocity ghost cells will be up-to-date
    EulerVel.FillBoundary(ParticleProperties::euler_velocity_index, 3, gm.periodicity());

    const int EulerVelocityIndex = ParticleProperties::euler_velocity_index;

    for(mParIter pti(*mContainer, LOCAL_LEVEL); pti.isValid(); ++pti){
        
        const Box& box = pti.validbox();
        
        auto& particles = pti.GetArrayOfStructs();
        auto *p_ptr = particles.data();
        const Long np = pti.numParticles();

        auto& attri = pti.GetAttribs();
        auto* Up = attri[P_ATTR::U_Marker].data();
        auto* Vp = attri[P_ATTR::V_Marker].data();
        auto* Wp = attri[P_ATTR::W_Marker].data();
        const auto& E = EulerVel.array(pti);

        amrex::ParallelFor(np, [=] 
        AMREX_GPU_DEVICE (int i) noexcept{
            VelocityInterpolation_cir(i, p_ptr[i], Up[i], Vp[i], Wp[i], E, EulerVelocityIndex, box.loVect(), box.hiVect(), plo, dx, type);
        });
    }
    if (verbose) mContainer->WriteAsciiFile(amrex::Concatenate("particle", 2));
    //amrex::Abort("stop here!");
}
#endif
#if (AMREX_SPACEDIM == 2)
void mFiber::VelocityInterpolation(MultiFab &EulerVel, DELTA_FUNCTION_TYPE type)
{
    if (verbose) amrex::Print() << "\tmFiber::VelocityInterpolation\n";

    const auto& gm = mContainer->GetParGDB()->Geom(LOCAL_LEVEL);
    auto plo = gm.ProbLoArray();
    auto dx = gm.CellSizeArray();

    // 边界填充
    EulerVel.FillBoundary(FiberProperties::euler_velocity_index, 2, gm.periodicity());


    const int EulerVelocityIndex = FiberProperties::euler_velocity_index;

    for (mFiberIter pti(*mContainer, LOCAL_LEVEL); pti.isValid(); ++pti) {
        const Box& box = pti.validbox();
        auto& particles = pti.GetArrayOfStructs();
        auto *p_ptr = particles.data();
        const Long np = pti.numParticles();

        auto& attri = pti.GetAttribs();
        auto* Up = attri[P_ATTR::U_Marker].data();
        auto* Vp = attri[P_ATTR::V_Marker].data();
        const auto& E = EulerVel.array(pti);

        amrex::ParallelFor(np, [=]
        AMREX_GPU_DEVICE (int i) noexcept {
            VelocityInterpolation_cir(i, p_ptr[i], Up[i], Vp[i], E, EulerVelocityIndex, box.loVect(), box.hiVect(), plo, dx, type);
        });
    }

    if (verbose) mContainer->WriteAsciiFile(amrex::Concatenate("fiber", 2));
}
#endif
#if (AMREX_SPACEDIM == 3)
template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ForceSpreading_cic (P const& p,
                         Real Px,
                         Real Py,
                         Real Pz,
                         ParticleReal& fxP,
                         ParticleReal& fyP,
                         ParticleReal& fzP,
                         ParticleReal& mxP,
                         ParticleReal& myP,
                         ParticleReal& mzP,
                         Array4<Real> const& E,
                         int EulerForceIndex,
                         Real dv,
                         GpuArray<Real,AMREX_SPACEDIM> const& plo,
                         GpuArray<Real,AMREX_SPACEDIM> const& dx,
                         DELTA_FUNCTION_TYPE type)
{
    //const Real d = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
    //plo to ii jj kk
    Real lx = (p.pos(0) - plo[0]) / dx[0];
    Real ly = (p.pos(1) - plo[1]) / dx[1];
    Real lz = (p.pos(2) - plo[2]) / dx[2];

    int i = static_cast<int>(Math::floor(lx));
    int j = static_cast<int>(Math::floor(ly));
    int k = static_cast<int>(Math::floor(lz));
    fxP *= dv;
    fyP *= dv;
    fzP *= dv;
    RealVect moment = RealVect((p.pos(0) - Px), (p.pos(1) - Py), (p.pos(2) - Pz)).crossProduct(RealVect(fxP, fyP, fzP));
    mxP = moment[0];
    myP = moment[1];
    mzP = moment[2];
    //lagrangian to Euler
    for(int ii = -2; ii < +3; ii++){
        for(int jj = -2; jj < +3; jj++){
            for(int kk = -2; kk < +3; kk ++){
                Real tU, tV, tW;
                const Real xi =plo[0] + (i + ii) * dx[0] + dx[0]/2;
                const Real yj =plo[1] + (j + jj) * dx[1] + dx[1]/2;
                const Real kz =plo[2] + (k + kk) * dx[2] + dx[2]/2;
                deltaFunction( p.pos(0), xi, dx[0], tU, type);
                deltaFunction( p.pos(1), yj, dx[1], tV, type);
                deltaFunction( p.pos(2), kz, dx[2], tW, type);
                Real delta_value = tU * tV * tW;
                Gpu::Atomic::AddNoRet(&E(i + ii, j + jj, k + kk, EulerForceIndex  ), delta_value * fxP);
                Gpu::Atomic::AddNoRet(&E(i + ii, j + jj, k + kk, EulerForceIndex+1), delta_value * fyP);
                Gpu::Atomic::AddNoRet(&E(i + ii, j + jj, k + kk, EulerForceIndex+2), delta_value * fzP);
            }
        }
    }
}
#endif
#if (AMREX_SPACEDIM == 2)
template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ForceSpreading_cic (P const& p,
                         Real Px, Real Py,
                         ParticleReal& fxP, ParticleReal& fyP,
                         ParticleReal& mxP, ParticleReal& myP,
                         Array4<Real> const& E,
                         int EulerForceIndex,
                         Real da,
                         GpuArray<Real,AMREX_SPACEDIM> const& plo,
                         GpuArray<Real,AMREX_SPACEDIM> const& dx,
                         DELTA_FUNCTION_TYPE type)
{
    // 坐标计算
    Real lx = (p.pos(0) - plo[0]) / dx[0];
    Real ly = (p.pos(1) - plo[1]) / dx[1];
    int i = static_cast<int>(Math::floor(lx));
    int j = static_cast<int>(Math::floor(ly));

    // 力乘以面积
    fxP *= da; fyP *= da;
    // 循环
    for(int ii = -2; ii < +3; ii++){
        for(int jj = -2; jj < +3; jj++){
            Real tU, tV;
            const Real xi = plo[0] + (i + ii) * dx[0] + dx[0]/2;
            const Real yj = plo[1] + (j + jj) * dx[1] + dx[1]/2;
            deltaFunction(p.pos(0), xi, dx[0], tU, type);
            deltaFunction(p.pos(1), yj, dx[1], tV, type);
            Real delta_value = tU * tV;
            Gpu::Atomic::AddNoRet(&E(i + ii, j + jj, 0, EulerForceIndex  ), delta_value * fxP);
            Gpu::Atomic::AddNoRet(&E(i + ii, j + jj, 0, EulerForceIndex+1), delta_value * fyP);
        }
    }
}
#endif

#if (AMREX_SPACEDIM == 3)
void mParticle::ForceSpreading(MultiFab & EulerForce,
                               kernel& kernel,
                               DELTA_FUNCTION_TYPE type)
{
    if (verbose) amrex::Print() << "\tmParticle::ForceSpreading\n";
    const auto& gm = mContainer->GetParGDB()->Geom(LOCAL_LEVEL);
    auto plo = gm.ProbLoArray();
    auto dxi = gm.CellSizeArray();
    for(mParIter pti(*mContainer, LOCAL_LEVEL); pti.isValid(); ++pti){
        const Long np = pti.numParticles();
        const auto& particles = pti.GetArrayOfStructs();
        auto Uarray = EulerForce[pti].array();
        auto& attri = pti.GetAttribs();

        auto *const fxP_ptr = attri[P_ATTR::Fx_Marker].data();
        auto *const fyP_ptr = attri[P_ATTR::Fy_Marker].data();
        auto *const fzP_ptr = attri[P_ATTR::Fz_Marker].data();
        auto *const mxP_ptr = attri[P_ATTR::Mx_Marker].data();
        auto *const myP_ptr = attri[P_ATTR::My_Marker].data();
        auto *const mzP_ptr = attri[P_ATTR::Mz_Marker].data();
        const auto *const p_ptr = particles().data();

        auto loc_ptr = kernel.location;
        auto dv = kernel.dv;
        auto force_index = ParticleProperties::euler_force_index;
        amrex::ParallelFor(np, [=]
        AMREX_GPU_DEVICE (int i) noexcept{
            ForceSpreading_cic(p_ptr[i], loc_ptr[0], loc_ptr[1], loc_ptr[2],
                               fxP_ptr[i], fyP_ptr[i], fzP_ptr[i], 
                               mxP_ptr[i], myP_ptr[i], mzP_ptr[i], 
                               Uarray, force_index, dv, plo, dxi, type);
        });
    }
    //barrier for sync;
    amrex::ParallelDescriptor::Barrier();

    using pc = mParticleContainer::SuperParticleType;
    // Each Processor
    auto fx = amrex::ReduceSum( *mContainer, [=]AMREX_GPU_HOST_DEVICE(const pc& p)->ParticleReal{return p.rdata(P_ATTR::Fx_Marker);});
    auto fy = amrex::ReduceSum( *mContainer, [=]AMREX_GPU_HOST_DEVICE(const pc& p)->ParticleReal{return p.rdata(P_ATTR::Fy_Marker);});
    auto fz = amrex::ReduceSum( *mContainer, [=]AMREX_GPU_HOST_DEVICE(const pc& p)->ParticleReal{return p.rdata(P_ATTR::Fz_Marker);});
    auto mx = amrex::ReduceSum( *mContainer, [=]AMREX_GPU_HOST_DEVICE(const pc& p)->ParticleReal{return p.rdata(P_ATTR::Mx_Marker);});
    auto my = amrex::ReduceSum( *mContainer, [=]AMREX_GPU_HOST_DEVICE(const pc& p)->ParticleReal{return p.rdata(P_ATTR::My_Marker);});
    auto mz = amrex::ReduceSum( *mContainer, [=]AMREX_GPU_HOST_DEVICE(const pc& p)->ParticleReal{return p.rdata(P_ATTR::Mz_Marker);});
    // MPI sum reduce -> current particle all IB force and moment
    amrex::ParallelAllReduce::Sum(fx, ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Sum(fy, ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Sum(fz, ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Sum(mx, ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Sum(my, ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Sum(mz, ParallelDescriptor::Communicator());

    kernel.ib_force = {fx, fy, fz};
    kernel.ib_moment = {mx, my, mz};

    EulerForce.SumBoundary(ParticleProperties::euler_force_index, 3, gm.periodicity());

    // if (false) {
    //     // Check the Multifab
    //     // Open a file stream for output
    //     std::ofstream outFile("EulerForce.txt");

    //     // Check the Multifab
    //     // for (MFIter mfi(EulerForce, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    //     for (MFIter mfi(EulerForce, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    //     {
    //         const Box& bx = mfi.validbox();
    //         outFile << "Box: " << bx << "\n"
    //                 << "From: (" << bx.smallEnd(0) << ", " << bx.smallEnd(1) << ", " << bx.smallEnd(2) << ") "
    //                 << "To: (" << bx.bigEnd(0) << ", " << bx.bigEnd(1) << ", " << bx.bigEnd(2) << ")\n";

    //         Array4<Real> const& a = EulerForce[mfi].array();

    //         // CPU context or illustrative purposes only
    //         for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
    //             for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
    //                 for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
    //                     // This print statement is for demonstration and should not be used in actual GPU code.
    //                     outFile << "Processing i: " << i << ", j: " << j << ", k: " << k << " " << a(i,j,k,0) << " " << a(i,j,k,1) << " " << a(i,j,k,2) << "\n";
    //                 }
    //             }
    //         }
    //     }

    //     // Close the file when done
    //     outFile.close();
    // }

}
#endif
#if (AMREX_SPACEDIM == 2)
void mFiber::ForceSpreading(MultiFab & EulerForce, fiber_kernel& kernel, DELTA_FUNCTION_TYPE type)
{
    if (verbose) amrex::Print() << "\tmFiber::ForceSpreading\n";

    const auto& gm = mContainer->GetParGDB()->Geom(LOCAL_LEVEL);
    auto plo = gm.ProbLoArray();
    auto dxi = gm.CellSizeArray();

    for(mFiberIter pti(*mContainer, LOCAL_LEVEL); pti.isValid(); ++pti){
        const Long np = pti.numParticles();
        const auto& particles = pti.GetArrayOfStructs();
        auto Uarray = EulerForce[pti].array();
        auto& attri = pti.GetAttribs();

        auto *const fxP_ptr = attri[P_ATTR::Fx_Marker].data();
        auto *const fyP_ptr = attri[P_ATTR::Fy_Marker].data();
        auto *const mxP_ptr = attri[P_ATTR::Mx_Marker].data();
        auto *const myP_ptr = attri[P_ATTR::My_Marker].data();
        const auto *const p_ptr = particles().data();

        // 细丝系统的中心坐标：使用起点作为参考点
        auto x0 = kernel.x_marker[0];  // 细丝起点 x 坐标
        auto y0 = kernel.y_marker[0];  // 细丝起点 y 坐标

        // 细丝系统的面积参数
        auto da = kernel.da;

        auto force_index = FiberProperties::euler_force_index;

        amrex::ParallelFor(np, [=]
        AMREX_GPU_DEVICE (int i) noexcept{
            ForceSpreading_cic(p_ptr[i], x0, y0,
                               fxP_ptr[i], fyP_ptr[i],
                               mxP_ptr[i], myP_ptr[i],
                               Uarray, force_index, da, plo, dxi, type);
        });
    }

    // 归约所有处理器的数据
    using pc = mFiberContainer::SuperParticleType;
    auto fx = amrex::ReduceSum(*mContainer, [=]AMREX_GPU_HOST_DEVICE(const pc& p)->ParticleReal{
        return p.rdata(P_ATTR::Fx_Marker);
    });
    auto fy = amrex::ReduceSum(*mContainer, [=]AMREX_GPU_HOST_DEVICE(const pc& p)->ParticleReal{
        return p.rdata(P_ATTR::Fy_Marker);
    });
    auto mx = amrex::ReduceSum(*mContainer, [=]AMREX_GPU_HOST_DEVICE(const pc& p)->ParticleReal{
        return p.rdata(P_ATTR::Mx_Marker);
    });
    auto my = amrex::ReduceSum(*mContainer, [=]AMREX_GPU_HOST_DEVICE(const pc& p)->ParticleReal{
        return p.rdata(P_ATTR::My_Marker);
    });

    // MPI归约
    amrex::ParallelAllReduce::Sum(fx, ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Sum(fy, ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Sum(mx, ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Sum(my, ParallelDescriptor::Communicator());

    // 二维下只有x和y方向的力
    kernel.ib_force = {fx, fy};

    // 边界同步
    EulerForce.SumBoundary(FiberProperties::euler_force_index, 2, gm.periodicity());
}
#endif
#if (AMREX_SPACEDIM == 3)
void mParticle::ResetLargrangianPoints()
{
    if (verbose) amrex::Print() << "\tmParticle::ResetLargrangianPoints\n";

    for(mParIter pti(*mContainer, LOCAL_LEVEL); pti.isValid(); ++pti){
        const Long np = pti.numParticles();
        auto& attri = pti.GetAttribs();

        auto *const vUP_ptr = attri[P_ATTR::U_Marker].data();
        auto *const vVP_ptr = attri[P_ATTR::V_Marker].data();
        auto *const vWP_ptr = attri[P_ATTR::W_Marker].data();
        auto *const fxP_ptr = attri[P_ATTR::Fx_Marker].data();
        auto *const fyP_ptr = attri[P_ATTR::Fy_Marker].data();
        auto *const fzP_ptr = attri[P_ATTR::Fz_Marker].data();
        auto *const mxP_ptr = attri[P_ATTR::Mx_Marker].data();
        auto *const myP_ptr = attri[P_ATTR::My_Marker].data();
        auto *const mzP_ptr = attri[P_ATTR::Mz_Marker].data();
        amrex::ParallelFor(np, [=]
        AMREX_GPU_DEVICE (int i) noexcept{
            vUP_ptr[i] = 0.0;
            vVP_ptr[i] = 0.0;
            vWP_ptr[i] = 0.0;
            fxP_ptr[i] = 0.0;
            fyP_ptr[i] = 0.0;
            fzP_ptr[i] = 0.0;
            mxP_ptr[i] = 0.0;
            myP_ptr[i] = 0.0;
            mzP_ptr[i] = 0.0;
        });
    }
}
#endif
#if (AMREX_SPACEDIM == 2)
void mFiber::ResetLargrangianPoints()
{
    if (verbose) amrex::Print() << "\tmFiber::ResetLargrangianPoints\n";

    for(mFiberIter pti(*mContainer, LOCAL_LEVEL); pti.isValid(); ++pti){
        const Long np = pti.numParticles();
        auto& attri = pti.GetAttribs();

        auto *const vUP_ptr = attri[P_ATTR::U_Marker].data();
        auto *const vVP_ptr = attri[P_ATTR::V_Marker].data();
        auto *const fxP_ptr = attri[P_ATTR::Fx_Marker].data();
        auto *const fyP_ptr = attri[P_ATTR::Fy_Marker].data();
        auto *const mxP_ptr = attri[P_ATTR::Mx_Marker].data();
        auto *const myP_ptr = attri[P_ATTR::My_Marker].data();
        
        amrex::ParallelFor(np, [=]
        AMREX_GPU_DEVICE (int i) noexcept{
            vUP_ptr[i] = 0.0;
            vVP_ptr[i] = 0.0;

            fxP_ptr[i] = 0.0;
            fyP_ptr[i] = 0.0;

            mxP_ptr[i] = 0.0;
            myP_ptr[i] = 0.0;

        });
    }
}
#endif
#if (AMREX_SPACEDIM == 3)
void mParticle::UpdateParticles(int iStep,
                                Real time,
                                const MultiFab& Euler_old, 
                                const MultiFab& Euler,
                                MultiFab& phi_nodal, 
                                MultiFab& pvf, 
                                Real dt)
{
    if (verbose) amrex::Print() << "mParticle::UpdateParticles\n";
    // start record
    auto UpdateParticlesStart = ParallelDescriptor::second();

    //Particle Collision calculation
    DoParticleCollision(ParticleProperties::collision_model);
    
    MultiFab AllParticlePVF(pvf.boxArray(), pvf.DistributionMap(), pvf.nComp(), pvf.nGrow());
    AllParticlePVF.setVal(0.0);
    
    //continue condition 6DOF
    for(auto& kernel : particle_kernels){

        calculate_phi_nodal(phi_nodal, kernel);
        nodal_phi_to_pvf(pvf, phi_nodal);

        // // fixed particle
        // if( ( kernel.TL.sum() == 0 ) &&
        //     ( kernel.RL.sum() == 0 ) ) {
        //     amrex::Print() << "Particle (" << kernel.id << ") is fixed\n";
        //     MultiFab::Add(AllParticlePVF, pvf, 0, 0, 1, 0); // do not copy ghost cell values
        //     continue;
        // }

        int ncomp = pvf.nComp();
        int ngrow = pvf.nGrow();
        MultiFab pvf_old(pvf.boxArray(), pvf.DistributionMap(), ncomp, ngrow);
        MultiFab::Copy(pvf_old, pvf, 0, 0, ncomp, ngrow);

        // bool at_least_one_free_trans_motion = ( kernel.TL[0] == 2 ) || 
        //                                       ( kernel.TL[1] == 2 ) ||
        //                                       ( kernel.TL[2] == 2 );
        // bool at_least_one_free_rot_motion   = ( kernel.RL[0] == 2 ) || 
        //                                       ( kernel.RL[1] == 2 ) ||
        //                                       ( kernel.RL[2] == 2 );

        int loop = ParticleProperties::loop_solid;

        while (loop > 0 && iStep > ParticleProperties::start_step) {

            // if(at_least_one_free_trans_motion) {
                kernel.sum_u_new.scale(0.0);
                kernel.sum_u_old.scale(0.0);
                // sum U
                CalculateSumU_cir(kernel.sum_u_new, Euler, pvf, ParticleProperties::euler_velocity_index);
                CalculateSumU_cir(kernel.sum_u_old, Euler_old, pvf_old, ParticleProperties::euler_velocity_index);
                amrex::ParallelAllReduce::Sum(kernel.sum_u_new.dataPtr(), 3, amrex::ParallelDescriptor::Communicator());
                amrex::ParallelAllReduce::Sum(kernel.sum_u_old.dataPtr(), 3, amrex::ParallelDescriptor::Communicator());
            // }

            // if(at_least_one_free_rot_motion) {
                kernel.sum_t_new.scale(0.0);
                kernel.sum_t_old.scale(0.0);
                // sum T
                CalculateSumT_cir(kernel.sum_t_new, Euler, pvf, kernel.location, ParticleProperties::euler_velocity_index);
                CalculateSumT_cir(kernel.sum_t_old, Euler_old, pvf_old, kernel.location, ParticleProperties::euler_velocity_index);
                amrex::ParallelAllReduce::Sum(kernel.sum_t_new.dataPtr(), 3, amrex::ParallelDescriptor::Communicator());
                amrex::ParallelAllReduce::Sum(kernel.sum_t_old.dataPtr(), 3, amrex::ParallelDescriptor::Communicator());
            // }

            // 6DOF
            if(ParallelDescriptor::MyProc() == ParallelDescriptor::IOProcessorNumber()){

                for(auto idir : {0,1,2})
                {
                    //TL
                    if (kernel.TL[idir] == 0) {
                        kernel.velocity[idir] = 0.0;
                    }
                    else if (kernel.TL[idir] == 1) {
                        kernel.location[idir] = kernel.location_old[idir] + (kernel.velocity[idir] + kernel.velocity_old[idir]) * dt * 0.5;
                    }
                    else if (kernel.TL[idir] == 2) {
                        if(!ParticleProperties::Uhlmann){
                            kernel.velocity[idir] = kernel.velocity_old[idir]
                                                + ((kernel.sum_u_new[idir] - kernel.sum_u_old[idir]) * ParticleProperties::euler_fluid_rho / dt 
                                                - kernel.ib_force[idir] * ParticleProperties::euler_fluid_rho
                                                + m_gravity[idir] * (kernel.rho - ParticleProperties::euler_fluid_rho) * kernel.Vp 
                                                + kernel.Fcp[idir]) * dt / kernel.rho / kernel.Vp ;
                        }else{
                            //Uhlmann
                            kernel.velocity[idir] = kernel.velocity_old[idir]
                                                + (ParticleProperties::euler_fluid_rho / kernel.Vp /(ParticleProperties::euler_fluid_rho - kernel.rho)*kernel.ib_force[idir]
                                                + m_gravity[idir]) * dt;
                        }
                        kernel.location[idir] = kernel.location_old[idir] + (kernel.velocity[idir] + kernel.velocity_old[idir]) * dt * 0.5;
                    }
                    else {
                        amrex::Print() << "Particle (" << kernel.id << ") has wrong TL"<< direction_str[idir] <<" value\n";
                        amrex::Abort("Stop here!");
                    }
                    //RL
                    if (kernel.RL[idir] == 0) {
                        kernel.omega[idir] = 0.0;
                    }
                    else if (kernel.RL[idir] == 1) {
                    }
                    else if (kernel.RL[idir] == 2) {
                        if(!ParticleProperties::Uhlmann){
                            kernel.omega[idir] = kernel.omega_old[idir]
                                            + ((kernel.sum_t_new[idir] - kernel.sum_t_old[idir]) * ParticleProperties::euler_fluid_rho / dt
                                            - kernel.ib_moment[idir] * ParticleProperties::euler_fluid_rho
                                            + kernel.Tcp[idir]) * dt / cal_momentum(kernel.rho, kernel.radius);
                        }else{
                            //Uhlmann
                            kernel.omega[idir] = kernel.omega_old[idir]
                                            + ParticleProperties::euler_fluid_rho /(ParticleProperties::euler_fluid_rho - kernel.rho) * kernel.ib_moment[idir] * kernel.dv
                                            / cal_momentum(kernel.rho, kernel.radius) * kernel.rho * dt;
                        }
                    }
                    else {
                        amrex::Print() << "Particle (" << kernel.id << ") has wrong RL"<< direction_str[idir] <<" value\n";
                        amrex::Abort("Stop here!");
                    }

                }
            }
            ParallelDescriptor::Bcast(&kernel.location[0],3,ParallelDescriptor::IOProcessorNumber());
            ParallelDescriptor::Bcast(&kernel.location_old[0],3,ParallelDescriptor::IOProcessorNumber());
            ParallelDescriptor::Bcast(&kernel.velocity[0],3,ParallelDescriptor::IOProcessorNumber());
            ParallelDescriptor::Bcast(&kernel.velocity_old[0],3,ParallelDescriptor::IOProcessorNumber());
            ParallelDescriptor::Bcast(&kernel.omega[0],3,ParallelDescriptor::IOProcessorNumber());
            ParallelDescriptor::Bcast(&kernel.omega_old[0],3,ParallelDescriptor::IOProcessorNumber());
        
            loop--;

            if (loop > 0) {
                calculate_phi_nodal(phi_nodal, kernel);
                nodal_phi_to_pvf(pvf, phi_nodal);
            }

        }
        
        RecordOldValue(kernel);
        MultiFab::Add(AllParticlePVF, pvf, 0, 0, 1, 0); // do not copy ghost cell values
    }
    // calculate the pvf based on the information of all particles
    MultiFab::Copy(pvf, AllParticlePVF, 0, 0, 1, pvf.nGrow());

    spend_time += ParallelDescriptor::second() - UpdateParticlesStart;
    ParallelDescriptor::ReduceRealMax(spend_time);
    amrex::Print() << "[DIBM] IB and update particle, step : "<< iStep <<", time : " << spend_time << "\n";

    int particle_write_freq = ParticleProperties::write_freq;
    if (iStep % particle_write_freq == 0) {
        for(auto kernel: particle_kernels) 
            WriteIBForceAndMoment(iStep, time, dt, kernel);
    }

    if (verbose) mContainer->WriteAsciiFile(amrex::Concatenate("particle", 4));
}
#endif
#if (AMREX_SPACEDIM == 2)
void mFiber::UpdateFibers(int iStep, Real time, Real dt)
{
    if (verbose) amrex::Print() << "mFiber::UpdateFibers\n";
            //time += dt;//时间推进
            for(auto& fiber_kernel : fiber_kernels){
            // 更新每个拉格朗日点的位置和速度
            for(int i = 0; i < fiber_kernel.num_marker; i++){
                Real x_pos = fiber_kernel.x_marker[i];
                Real y_pos = fiber_kernel.y0_base + fiber_kernel.amplitude * 
                            sin(2 * M_PI * (time / fiber_kernel.period - fiber_kernel.wave_number * x_pos / fiber_kernel.length) + fiber_kernel.phase);
                
                // 更新位置
                fiber_kernel.y_marker[i] = y_pos;
                
                // 计算目标速度
                Real u_b_y = (2 * M_PI * fiber_kernel.amplitude / fiber_kernel.period) * 
                            cos(2 * M_PI * (time / fiber_kernel.period - fiber_kernel.wave_number * x_pos / fiber_kernel.length) + fiber_kernel.phase);
                fiber_kernel.velocity_marker[i] = RealVect(0.0, u_b_y);
            }
            
            fiber_kernel.y_marker_d.resize(fiber_kernel.y_marker.size());
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                                  fiber_kernel.y_marker.begin(), fiber_kernel.y_marker.end(),
                                  fiber_kernel.y_marker_d.begin());

            fiber_kernel.velocity_marker_d.resize(fiber_kernel.velocity_marker.size());
            amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice,
                                  fiber_kernel.velocity_marker.begin(), fiber_kernel.velocity_marker.end(),
                                  fiber_kernel.velocity_marker_d.begin());
            amrex::Gpu::streamSynchronize();
        }

        // 输出拉格朗日点到 CSV 文件,可视化细丝变化
        //std::vector<std::string> all_data;
    //for (const auto& fiber_kernel : fiber_kernels) {
    //    std::string file_name = "Fiber_LagrangianPoints_" + std::to_string(fiber_kernel.id) + ".csv";
    //    std::ofstream out_file;

        // 如果文件不存在，写入表头
    //    if (!fs::exists(file_name)) {
    //        out_file.open(file_name, std::ios::app);
    //        out_file << "Step,Time";
    //        for (int i = 0; i < fiber_kernel.num_marker; i++) {
    //            out_file << ",x" << i << ",y" << i;
    //        }
    //        out_file << "\n";
    //    } else {
    //        out_file.open(file_name, std::ios::app);
    //    }

        // 写入当前步的数据
     //   out_file << iStep << "," << time;
     //   for (int i = 0; i < fiber_kernel.num_marker; i++) {
     //       out_file << "," << fiber_kernel.x_marker[i] << "," << fiber_kernel.y_marker[i];
     //   }
     //   out_file << "\n";
     //   out_file.close();
    //}

    int fiber_write_freq = FiberProperties::write_freq;
    if (iStep % fiber_write_freq == 0) {
        for(auto kernel: fiber_kernels) 
            WriteIBForceAndMoment(iStep, time, dt, kernel);
    }
    if (verbose) mContainer->WriteAsciiFile(amrex::Concatenate("fiber", 4));
}
#endif
#if (AMREX_SPACEDIM == 3)
void mParticle::DoParticleCollision(int model)
{
    if(particle_kernels.size() < 2 ) return ;

    if (verbose) amrex::Print() << "\tmParticle::DoParticleCollision\n";
    
    if(ParallelDescriptor::MyProc() == ParallelDescriptor::IOProcessorNumber()){
        for(const auto& kernel : particle_kernels){
            m_Collision.InsertParticle(kernel.location, kernel.velocity, kernel.radius, kernel.rho);
        }
        
        m_Collision.takeModel(model);

        for(auto & particle_kernel : particle_kernels){
            particle_kernel.Fcp = m_Collision.Particles.front().preForece 
                                * particle_kernel.Vp * particle_kernel.rho * m_gravity.vectorLength();
            m_Collision.Particles.pop_front();
        }
    }
    for(auto& kernel : particle_kernels){
        ParallelDescriptor::Bcast(kernel.Fcp.dataPtr(), 3, ParallelDescriptor::IOProcessorNumber());
    }
}
#endif
#if (AMREX_SPACEDIM == 3)
void mParticle::ComputeLagrangianForce(Real dt, 
                                       const kernel& kernel)
{
    
    if (verbose) amrex::Print() << "\tmParticle::ComputeLagrangianForce\n";

    Real Ub = kernel.velocity[0];
    Real Vb = kernel.velocity[1];
    Real Wb = kernel.velocity[2];
    Real Px = kernel.location[0];
    Real Py = kernel.location[1];
    Real Pz = kernel.location[2];

    for(mParIter pti(*mContainer, LOCAL_LEVEL); pti.isValid(); ++pti){
        const Long np = pti.numParticles();
        auto& attri = pti.GetAttribs();
        auto const* p_ptr = pti.GetArrayOfStructs().data();

        auto* Up = attri[P_ATTR::U_Marker].data();
        auto* Vp = attri[P_ATTR::V_Marker].data();
        auto* Wp = attri[P_ATTR::W_Marker].data();
        auto *FxP = attri[P_ATTR::Fx_Marker].data();
        auto *FyP = attri[P_ATTR::Fy_Marker].data();
        auto *FzP = attri[P_ATTR::Fz_Marker].data();

        amrex::ParallelFor(np,
        [=] AMREX_GPU_DEVICE (int i) noexcept{
            auto Ur = (kernel.omega).crossProduct(RealVect(p_ptr[i].pos(0) - Px, p_ptr[i].pos(1) - Py, p_ptr[i].pos(2) - Pz));
            FxP[i] = (Ub + Ur[0] - Up[i])/dt; //
            FyP[i] = (Vb + Ur[1] - Vp[i])/dt; //
            FzP[i] = (Wb + Ur[2] - Wp[i])/dt; //
        });
    }
    if (verbose) mContainer->WriteAsciiFile(amrex::Concatenate("particle", 3));
}
#endif
#if (AMREX_SPACEDIM == 2)
void mFiber::ComputeLagrangianForce(Real dt, const fiber_kernel& kernel)
{
    if (verbose) amrex::Print() << "\tmFiber::ComputeLagrangianForce\n";
    
    for(mFiberIter pti(*mContainer, LOCAL_LEVEL); pti.isValid(); ++pti){
        const Long np = pti.numParticles();
        auto& attri = pti.GetAttribs();
        auto const* p_ptr = pti.GetArrayOfStructs().data();
        
        auto* Up = attri[P_ATTR::U_Marker].data();
        auto* Vp = attri[P_ATTR::V_Marker].data();
        auto *FxP = attri[P_ATTR::Fx_Marker].data();
        auto *FyP = attri[P_ATTR::Fy_Marker].data();
        //const RealVect* velocity_marker_ptr = kernel.velocity_marker.data();
        const RealVect* velocity_marker_ptr = kernel.velocity_marker_d.data();

        // 系数（可改成输入参数）
        const Real k_pos = 1000.0;
        const Real k_vel = 1.0;
        //const Real k_penalty = 1000.0;

        //const Real* xex = kernel.x_marker.data();
        //const Real* xey = kernel.y_marker.data();
        const Real* xex = kernel.x_marker_d.data();
        const Real* xey = kernel.y_marker_d.data();
        amrex::ParallelFor(np,
        [=] AMREX_GPU_DEVICE (int i) noexcept{
// IB力由速度差计算*********************************************************************************
            //const int id = p_ptr[i].id() - 1; // 粒子 id 从 1 开始，这里转换为 0-based 索引
            //Real u_b_x = 0.0;
            //Real u_b_y = fiber_kernel.velocity_marker[i][1]; // 从预计算的速度获取
            //Real u_b_y = velocity_marker_ptr[id][1]; // 从预计算的速度获取
            //计算IB力
            //FxP[i] = (u_b_x - Up[i]) / dt;
            //FyP[i] = (u_b_y - Vp[i]) / dt;

//由位置差计算IB力*********************************************************************************
            // X: 当前拉格朗日点位置
        //    const Real Xx = p_ptr[i].pos(0);
        //    const Real Xy = p_ptr[i].pos(1);

            // id 从 1 开始
        //    const int id = p_ptr[i].id() - 1;

            // Xe: 期望位置（由运动公式生成并保存在 kernel.x_marker/y_marker）
        //    const Real Xex = kernel.x_marker[id];
        //    const Real Xey = kernel.y_marker[id];

            // F = k (Xe - X)
        //    FxP[i] = k_penalty * (Xex - Xx);
        //    FyP[i] = k_penalty * (Xey - Xy);
//IB力由位置差和速度差共同决定***********************************************************************
            const int id = p_ptr[i].id() - 1; // 粒子 id 从 1 开始，这里转换为 0-based 索引

            // X（当前真实位置）与 Xe（期望位置）
            const Real Xx  = p_ptr[i].pos(0);
            const Real Xy  = p_ptr[i].pos(1);
            const Real Xex = xex[id];
            const Real Xey = xey[id];

            // U（插值得到的真实速度）与 Ue（期望速度）
            const Real Ux  = Up[i];
            const Real Uy  = Vp[i];
            const Real Uex = 0.0;
            const Real Uey = velocity_marker_ptr[id][1];

            // 合成惩罚力：F = k_pos (Xe - X) + k_vel (Ue - U) / dt
            FxP[i] = k_pos * (Xex - Xx) / dt + k_vel * (Uex - Ux) / dt;
            FyP[i] = k_pos * (Xey - Xy) / dt + k_vel * (Uey - Uy) / dt;
        });
    }
    if (verbose) mContainer->WriteAsciiFile(amrex::Concatenate("fiber", 3));
}
#endif
#if (AMREX_SPACEDIM == 3)
void mParticle::VelocityCorrection(amrex::MultiFab &Euler, amrex::MultiFab &EulerForce, Real dt) const
{
    if(verbose) amrex::Print() << "\tmParticle::VelocityCorrection\n";
    MultiFab::Saxpy(Euler, dt, EulerForce, ParticleProperties::euler_force_index, ParticleProperties::euler_velocity_index, 3, 0); //VelocityCorrection
}
#endif
#if (AMREX_SPACEDIM == 2)
void mFiber::VelocityCorrection(amrex::MultiFab &Euler, amrex::MultiFab &EulerForce, Real dt) const
{
    if(verbose) amrex::Print() << "\tmFiber::VelocityCorrection\n";
    
    // 速度修正
    MultiFab::Saxpy(Euler, dt, EulerForce, FiberProperties::euler_force_index, FiberProperties::euler_velocity_index, 2, 0);

}
#endif
#if (AMREX_SPACEDIM == 3)
void mParticle::RecordOldValue(kernel& kernel)
{
    kernel.location_old = kernel.location;
    kernel.velocity_old = kernel.velocity;
    kernel.omega_old = kernel.omega;
}
#endif
#if (AMREX_SPACEDIM == 3)
void mParticle::WriteParticleFile(int index)
{
    mContainer->WriteAsciiFile(amrex::Concatenate("particle", index));
}
#endif
#if (AMREX_SPACEDIM == 2)
void mFiber::WriteFiberFile(int index)
{
    mContainer->WriteAsciiFile(amrex::Concatenate("fiber", index));
}
#endif
#if (AMREX_SPACEDIM == 3)
void mParticle::WriteIBForceAndMoment(int step, amrex::Real time, amrex::Real dt, kernel& current_kernel)
{
    
    if(amrex::ParallelDescriptor::MyProc() != ParallelDescriptor::IOProcessorNumber()) return; 

    std::string file("IB_Particle_" + std::to_string(current_kernel.id) + ".csv");
    std::ofstream out_ib_force;

    std::string head;
    if(!fs::exists(file)){
        head = "iStep,time,X,Y,Z,Vx,Vy,Vz,Rx,Ry,Rz,Fx,Fy,Fz,Mx,My,Mz,Fcpx,Fcpy,Fcpz,Tcpx,Tcpy,Tcpz,SumUx,SumUy,SumUz,SumTx,SumTy,SumTz\n";
    }else{
        head = "";
    }

    out_ib_force.open(file, std::ios::app);
    if(!out_ib_force.is_open()){
        amrex::Print() << "[Particle] write particle file error , step: " << step;
    }else{
        out_ib_force << head << step << "," << time << "," 
                     << current_kernel.location[0] << "," << current_kernel.location[1] << "," << current_kernel.location[2] << ","
                     << current_kernel.velocity[0] << "," << current_kernel.velocity[1] << "," << current_kernel.velocity[2] << ","
                     << current_kernel.omega[0] << "," << current_kernel.omega[1] << "," << current_kernel.omega[2] << ","
                     << current_kernel.ib_force[0] << "," << current_kernel.ib_force[1] << "," << current_kernel.ib_force[2] << "," 
                     << current_kernel.ib_moment[0] << "," << current_kernel.ib_moment[1] << "," << current_kernel.ib_moment[2] << ","
                     << current_kernel.Fcp[0] << "," << current_kernel.Fcp[1] << "," << current_kernel.Fcp[2] << ","
                     << current_kernel.Tcp[0] << "," << current_kernel.Tcp[1] << "," << current_kernel.Tcp[2] << ","
                     << (current_kernel.sum_u_new[0] - current_kernel.sum_u_old[0])/dt << ","
                     << (current_kernel.sum_u_new[1] - current_kernel.sum_u_old[1])/dt << ","
                     << (current_kernel.sum_u_new[2] - current_kernel.sum_u_old[2])/dt << ","
                     << (current_kernel.sum_t_new[0] - current_kernel.sum_t_old[0])/dt << ","
                     << (current_kernel.sum_t_new[1] - current_kernel.sum_t_old[1])/dt << ","
                     << (current_kernel.sum_t_new[2] - current_kernel.sum_t_old[2])/dt << "\n";
    }
    out_ib_force.close();
}
#endif
#if (AMREX_SPACEDIM == 2)
void mFiber::WriteIBForceAndMoment(int step, amrex::Real time, amrex::Real dt, fiber_kernel& current_kernel)
{
    if(amrex::ParallelDescriptor::MyProc() != ParallelDescriptor::IOProcessorNumber()) return;

    std::string file("IB_Fiber_" + std::to_string(current_kernel.id) + ".csv");
    std::ofstream out_ib_force;

    std::string head;
    if(!fs::exists(file)){
        head = "iStep,time,Fx,Fy\n";
    }else{
        head = "";
    }

    out_ib_force.open(file, std::ios::app);
    if(!out_ib_force.is_open()){
        amrex::Print() << "[Fiber] write fiber file error , step: " << step;
    }else{
        out_ib_force << head << step << "," << time << ","
                     << current_kernel.ib_force[0] << "," << current_kernel.ib_force[1] << "\n";
    }
    out_ib_force.close();
}
#endif
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                    Particles member function                  */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#if (AMREX_SPACEDIM == 3)
void Particles::create_particles(const Geometry &gm,
                                 const DistributionMapping & dm,
                                 const BoxArray & ba)
{
    amrex::Print() << "[Particle] : create Particle Container\n";
    if(particle->mContainer != nullptr){
        delete particle->mContainer;
        particle->mContainer = nullptr;
    }
    particle->mContainer = new mParticleContainer(gm, dm, ba);

    //get particle tile
    std::pair<int, int> key{0,0};
    auto& particleTileTmp = particle->mContainer->GetParticles(0)[key];
    //insert markers
    if ( ParallelDescriptor::MyProc() == ParallelDescriptor::IOProcessorNumber() ) {
        //insert particle's markers
        //Real phiK = 0;
        for(int marker_index = 0; marker_index < particle->particle_kernels[0].ml; marker_index++){
            //insert code
            mParticleContainer::ParticleType markerP;
            markerP.id() = marker_index + 1;
            markerP.cpu() = ParallelDescriptor::MyProc();
            markerP.pos(0) = particle->particle_kernels[0].location[0];
            markerP.pos(1) = particle->particle_kernels[0].location[1];
            markerP.pos(2) = particle->particle_kernels[0].location[2];

            std::array<ParticleReal, numAttri> Marker_attr;
            Marker_attr[U_Marker] = 0.0;
            Marker_attr[V_Marker] = 0.0;
            Marker_attr[W_Marker] = 0.0;
            Marker_attr[Fx_Marker] = 0.0;
            Marker_attr[Fy_Marker] = 0.0;
            Marker_attr[Fz_Marker] = 0.0;

            particleTileTmp.push_back(markerP);
            particleTileTmp.push_back_real(Marker_attr);
        }
    }
    particle->mContainer->Redistribute(); // Still needs to redistribute here! 

    ParticleProperties::plo = gm.ProbLoArray();
    ParticleProperties::phi = gm.ProbHiArray();
    ParticleProperties::dx = gm.CellSizeArray();
}

mParticle* Particles::get_particles()
{
    return particle;
}
#endif
#if (AMREX_SPACEDIM == 2)
void Fibers::create_fibers(const Geometry &gm,
                           const DistributionMapping &dm,
                           const BoxArray &ba)
{
    amrex::Print() << "[Fiber] : create Fiber Container\n";
    if (fiber->mContainer != nullptr) {
        delete fiber->mContainer;
        fiber->mContainer = nullptr;
    }
    fiber->mContainer = new mFiberContainer(gm, dm, ba);

    // 获取fiber tile
    std::pair<int, int> key{0, 0};
    auto& fiberTileTmp = fiber->mContainer->GetParticles(0)[key];

    // 插入fiber节点
    if (ParallelDescriptor::MyProc() == ParallelDescriptor::IOProcessorNumber()) {
        for (int fiber_index = 0; fiber_index < fiber->fiber_kernels[0].num_marker; ++fiber_index) {
            mFiberContainer::ParticleType fiberNode;
            fiberNode.id() = fiber_index + 1;
            fiberNode.cpu() = ParallelDescriptor::MyProc();
            fiberNode.pos(0) = fiber->fiber_kernels[0].x_marker[fiber_index];
            fiberNode.pos(1) = fiber->fiber_kernels[0].y_marker[fiber_index];
            //fiberNode.pos(2) = 0.0; // 二维细丝，z坐标为0

            std::array<ParticleReal, numAttri> Fiber_attr;
            Fiber_attr[U_Marker] = fiber->fiber_kernels[0].velocity_marker[fiber_index][0];
            Fiber_attr[V_Marker] = fiber->fiber_kernels[0].velocity_marker[fiber_index][1];
            Fiber_attr[W_Marker] = 0.0; // 二维细丝，W分量为0
            Fiber_attr[Fx_Marker] = 0.0;
            Fiber_attr[Fy_Marker] = 0.0;
            Fiber_attr[Fz_Marker] = 0.0; // 二维细丝，Fz分量为0

            fiberTileTmp.push_back(fiberNode);
            fiberTileTmp.push_back_real(Fiber_attr);
        }
    }
    fiber->mContainer->Redistribute();

    //FiberProperties::plo = gm.ProbLoArray();
    auto prob_lo = gm.ProbLoArray();
    FiberProperties::plo = amrex::GpuArray<double,2>{prob_lo[0], prob_lo[1]};
    //FiberProperties::phi = gm.ProbHiArray();
    auto prob_hi = gm.ProbHiArray();
    FiberProperties::phi = amrex::GpuArray<double,2>{prob_hi[0], prob_hi[1]};
    //FiberProperties::dx = gm.CellSizeArray();
    auto cell_size = gm.CellSizeArray();
    FiberProperties::dx = amrex::GpuArray<double,2>{cell_size[0], cell_size[1]};
}
mFiber* Fibers::get_fibers()
{
    return fiber;
}
#endif
#if (AMREX_SPACEDIM == 3)
void Particles::init_particle(Real gravity, Real h)
{
    amrex::Print() << "[Particle] : create Particle's kernel\n";
    particle = new mParticle;
    if(particle != nullptr){
        isInitial = true;
        particle->InitParticles(
            ParticleProperties::_x, 
            ParticleProperties::_y, 
            ParticleProperties::_z,
            ParticleProperties::_rho,
            ParticleProperties::Vx,
            ParticleProperties::Vy,
            ParticleProperties::Vz,
            ParticleProperties::Ox,
            ParticleProperties::Oy,
            ParticleProperties::Oz,
            ParticleProperties::TLX,
            ParticleProperties::TLY,
            ParticleProperties::TLZ,
            ParticleProperties::RLX,
            ParticleProperties::RLY,
            ParticleProperties::RLZ,
            ParticleProperties::_radius,
            h,
            gravity,
            ParticleProperties::verbose);
    }

}
#endif
#if (AMREX_SPACEDIM == 2)
void Fibers::init_fiber(Real gravity, Real h)
{
    amrex::Print() << "[Fiber] : create Fiber's kernel\n";
    fiber = new mFiber;
    if (fiber != nullptr) {
        isInitial = true;
        fiber->InitFibers(
            FiberProperties::x0,
            FiberProperties::y0,
            FiberProperties::amplitude,
            FiberProperties::period,
            FiberProperties::length,
            FiberProperties::wave_number,
            FiberProperties::phase,
            FiberProperties::num_marker,
            h,
            FiberProperties::verbose);
    }
}
#endif
#if (AMREX_SPACEDIM == 3)
void Particles::Restart(Real gravity, Real h, int iStep)
{
    amrex::Print() << "[Particle] : restart Particle's kernel, step :" << iStep << "\n"
                   << "\tstart read particle csv file , default name is IB_Particle_x.csv\n" 
                   << "\tdo not delete those file before \"restart\"\n\n";
    delete particle;
    particle = new mParticle;
            particle->InitParticles(
            ParticleProperties::_x, 
            ParticleProperties::_y, 
            ParticleProperties::_z,
            ParticleProperties::_rho,
            ParticleProperties::Vx,
            ParticleProperties::Vy,
            ParticleProperties::Vz,
            ParticleProperties::Ox,
            ParticleProperties::Oy,
            ParticleProperties::Oz,
            ParticleProperties::TLX,
            ParticleProperties::TLY,
            ParticleProperties::TLZ,
            ParticleProperties::RLX,
            ParticleProperties::RLY,
            ParticleProperties::RLZ,
            ParticleProperties::_radius,
            h,
            gravity,
            ParticleProperties::verbose);
    //deal in IO processor
    //start read csv file
    for(auto& kernel : particle->particle_kernels){
        //filename
        if(amrex::ParallelDescriptor::MyProc() == amrex::ParallelDescriptor::IOProcessorNumber()){
            std::string fileName = "IB_Particle_" + std::to_string(kernel.id) + ".csv";
            std::string tmpfile = "tmp" + fileName;
            //file stream
            std::ifstream particle_data(fileName);
            std::ofstream particle_file(tmpfile);
            // open state
            if(!particle_data.is_open() || !particle_file.is_open()){
                amrex::Abort("\tCan not open particle file : " + fileName);
            }
            std::string lineData;
            int line{0};
            while(std::getline(particle_data, lineData)){
                line++;
                particle_file << lineData << "\n";
                if(line <= iStep) {
                    continue;
                }
                //old location
                //iStep,time,X,Y,Z,Vx,Vy,Vz,Rx,Ry,Rz,Fx,Fy,Fz,Mx,My,Mz,Fcpx,Fcpy,Fcpz,Tcpx,Tcpy,Tcpz
                if(line == iStep + 1){
                    std::stringstream ss(lineData);
                    std::string data;
                    std::vector<amrex::Real> dataStruct;
                    while(std::getline(ss, data, ',')){
                        dataStruct.emplace_back(std::stod(data));
                    }
                    kernel.location[0] = dataStruct[2];
                    kernel.location[1] = dataStruct[3];
                    kernel.location[2] = dataStruct[4];
                    kernel.velocity[0] = dataStruct[5];
                    kernel.velocity[1] = dataStruct[6];
                    kernel.velocity[2] = dataStruct[7];
                    kernel.omega[0] = dataStruct[8];
                    kernel.omega[1] = dataStruct[9];
                    kernel.omega[2] = dataStruct[10];

                    kernel.location_old = kernel.location;
                    kernel.velocity_old = kernel.velocity;
                    kernel.omega_old    = kernel.omega;
                    break;
                }
                else
                    break;
            }
            particle_data.close();
            particle_file.close();
            std::remove(fileName.c_str());
            std::rename(tmpfile.c_str(), fileName.c_str());
        }
        ParallelDescriptor::Bcast(&kernel.location[0], 3, ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::Bcast(&kernel.velocity[0], 3,ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::Bcast(&kernel.omega[0], 3,ParallelDescriptor::IOProcessorNumber());

        ParallelDescriptor::Bcast(&kernel.location_old[0], 3, ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::Bcast(&kernel.velocity_old[0], 3,ParallelDescriptor::IOProcessorNumber());
        ParallelDescriptor::Bcast(&kernel.omega_old[0], 3,ParallelDescriptor::IOProcessorNumber());
    }

    isInitial = true;
}
#endif
#if (AMREX_SPACEDIM == 3)
void Particles::Initialize()
{
    ParmParse pp("particle");

    std::string particle_inputfile;
    std::string particle_init_file;
    pp.get("input",particle_inputfile);
    
    if(!particle_inputfile.empty()){
        ParmParse p_file(particle_inputfile);
        p_file.query("init", particle_init_file);
        p_file.getarr("x",          ParticleProperties::_x);
        p_file.getarr("y",          ParticleProperties::_y);
        p_file.getarr("z",          ParticleProperties::_z);
        p_file.getarr("rho",        ParticleProperties::_rho);
        p_file.getarr("velocity_x", ParticleProperties::Vx);
        p_file.getarr("velocity_y", ParticleProperties::Vy);
        p_file.getarr("velocity_z", ParticleProperties::Vz);
        p_file.getarr("omega_x",    ParticleProperties::Ox);
        p_file.getarr("omega_y",    ParticleProperties::Oy);
        p_file.getarr("omega_z",    ParticleProperties::Oz);
        p_file.getarr("TLX",        ParticleProperties::TLX);
        p_file.getarr("TLY",        ParticleProperties::TLY);
        p_file.getarr("TLZ",        ParticleProperties::TLZ);
        p_file.getarr("RLX",        ParticleProperties::RLX);
        p_file.getarr("RLY",        ParticleProperties::RLY);
        p_file.getarr("RLZ",        ParticleProperties::RLZ);
        p_file.getarr("radius",     ParticleProperties::_radius);
        p_file.query("RD",          ParticleProperties::rd);
        p_file.query("LOOP_NS",     ParticleProperties::loop_ns);
        p_file.query("LOOP_SOLID",  ParticleProperties::loop_solid);
        p_file.query("verbose",     ParticleProperties::verbose);
        p_file.query("start_step",  ParticleProperties::start_step);
        p_file.query("Uhlmann",     ParticleProperties::Uhlmann);
        p_file.query("collision_model", ParticleProperties::collision_model);
        p_file.query("write_freq",  ParticleProperties::write_freq);
        
        ParmParse ns("ns");
        ns.get("fluid_rho",      ParticleProperties::euler_fluid_rho);
        
        ParmParse level_parse("amr");
        level_parse.get("max_level", ParticleProperties::euler_finest_level);

        ParmParse geometry_parse("geometry");
        geometry_parse.getarr("prob_lo", ParticleProperties::GLO);
        geometry_parse.getarr("prob_hi", ParticleProperties::GHI);
        amrex::Print() << "[Particle] : Reading partilces cfg file : " << particle_inputfile << "\n"
                       << "             Particle's level : " << ParticleProperties::euler_finest_level << "\n";

        if(!particle_init_file.empty()){
            ParticleProperties::init_particle_from_file = true;
            //clear particle position container
            ParticleProperties::_x.clear();
            ParticleProperties::_y.clear();
            ParticleProperties::_z.clear();
            // parse particle's location data
            std::ifstream init_particle(particle_init_file);
            std::string line_data;
            while(std::getline(init_particle, line_data)){
                // id x_location y_location z_location
                std::istringstream line(line_data);
                std::vector<std::string> str_tokne;
                std::string token;
                while(line >> token){
                    str_tokne.push_back(token);
                }

                ParticleProperties::_x.push_back(std::stod(str_tokne[0]));
                ParticleProperties::_y.push_back(std::stod(str_tokne[1]));
                ParticleProperties::_z.push_back(std::stod(str_tokne[2]));
            }
            ParticleProperties::_x.shrink_to_fit();
            ParticleProperties::_y.shrink_to_fit();
            ParticleProperties::_z.shrink_to_fit();
            amrex::Print() << "             initial Particle by file : " << particle_init_file
                           << "             particle's size : " << ParticleProperties::_x.size() << "\n";
        }

    }else {
        amrex::Abort("[Particle] : can't read particles settings, pls check your config file \"particle.input\"");
    }
}
#endif
#if (AMREX_SPACEDIM == 2)
void Fibers::Initialize()
{
    ParmParse pp("fiber");

    std::string fiber_inputfile;
    pp.get("input", fiber_inputfile);
    
    if(!fiber_inputfile.empty()){
        ParmParse p_file(fiber_inputfile);
        p_file.getarr("x0",          FiberProperties::x0);
        p_file.getarr("y0",          FiberProperties::y0);
        p_file.getarr("amplitude",   FiberProperties::amplitude);
        p_file.getarr("period",      FiberProperties::period);
        p_file.getarr("length",      FiberProperties::length);
        p_file.getarr("wave_number", FiberProperties::wave_number);
        p_file.getarr("phase",       FiberProperties::phase);
        p_file.getarr("num_marker",  FiberProperties::num_marker);
        p_file.query("verbose",      FiberProperties::verbose);
        p_file.query("start_step",   FiberProperties::start_step);
        p_file.query("write_freq",   FiberProperties::write_freq);
        p_file.query("LOOP_NS",      FiberProperties::loop_ns);
        
        ParmParse ns("ns");
        ns.get("fluid_rho",          FiberProperties::euler_fluid_rho);
        
        ParmParse level_parse("amr");
        level_parse.get("max_level", FiberProperties::euler_finest_level);

        ParmParse geometry_parse("geometry");
        geometry_parse.getarr("prob_lo", FiberProperties::GLO);
        geometry_parse.getarr("prob_hi", FiberProperties::GHI);
        
        amrex::Print() << "[Fiber] : Reading fiber cfg file : " << fiber_inputfile << "\n"
                       << "             Fiber's level : " << FiberProperties::euler_finest_level << "\n";

    }else {
        amrex::Abort("[Fiber] : can't read fiber settings, pls check your config file \"fiber.input\"");
    }
}
#endif
#if (AMREX_SPACEDIM == 3)
int Particles::ParticleFinestLevel()
{
    return ParticleProperties::euler_finest_level;
}
#endif
#if (AMREX_SPACEDIM == 2)
int Fibers::FiberFinestLevel()
{
    return FiberProperties::euler_finest_level;
}
#endif