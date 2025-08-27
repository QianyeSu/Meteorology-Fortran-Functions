# Fortran代码功能分析

## 1. KolmSmir2.f [📝 src](fortran/KolmSmir2.f)
**功能**: 实现Kolmogorov-Smirnov双样本检验

**详细说明**:
- 用于检验两个样本分布函数的差异
- 包含两个子程序:
  - KOLM2: 主检验程序，对输入的两个样本X和Y进行统计检验
  - SMIRN: 计算统计量对应的概率值
- 输入参数:
  - X: 样本1数据向量(长度N)
  - Y: 样本2数据向量(长度M) 
  - N,M: 两个样本的大小
  - IFLAG: 控制标志
- 输出参数:
  - Z: 检验统计量
  - PROB: 对应的概率值
  - D: 两个分布函数的最大差值
- 功能特点:
  - 自动对输入数据进行升序排序
  - 计算经验分布函数之间的最大差异
  - 适用于大样本检验(N,M>=100)
  - 基于经典的Kolmogorov-Smirnov统计理论

## 2. VSPCORR.f [📝 src](fortran/VSPCORR.f)  
**功能**: 计算Spearman等级相关系数(处理缺失值版本)

**详细说明**:
- 包含四个子程序:
  - SPCORRZ: 主入口程序，处理包含缺失值的数据
  - SPCORR: 计算Spearman等级相关系数的核心程序
  - SPRANK: 对数据进行等级排序
  - SPSORT: 二进制排序算法
- 输入参数:
  - X,Y: 输入数据向量(长度N)
  - N: 样本数量
  - HASXMSG,HASYMSG: 缺失值标志
  - XMSG,YMSG: 缺失值代码
  - IWRITE: 输出控制参数
- 输出参数:
  - SPC: Spearman等级相关系数(-1到1之间)
  - NMSG: 缺失值数量
- 功能特点:
  - 自动处理和排除缺失值
  - 使用高效的二进制排序算法
  - 支持相同值的平均等级处理
  - 适用于非参数相关分析
  - 基于经典的Spearman等级相关理论

## 3. aam.f [📝 src](fortran/aam.f)
**功能**: 计算大气相对角动量(Atmospheric Angular Momentum)

**详细说明**:
- 包含三个子程序:
  - DAAMOM1: 计算大气角动量(使用1D气压厚度数组)
  - DAAMOM3: 计算大气角动量(使用3D气压厚度数组) 
  - DCOSWEIGHT: 计算纬度权重
- 输入参数:
  - U: 纬向风速分量(m/s)，3维数组(经度×纬度×层数)
  - DP: 气压厚度(Pa)，可为1D或3D数组
  - LAT: 纬度数组
  - WGT: 纬度权重数组
  - MLON,NLAT,KLVL: 经度、纬度、垂直层数
  - UMSG: 风速缺失值代码
- 输出参数:
  - AAM: 大气相对角动量(kg$\cdot$m$^2$/s)
- 功能特点:
  - 基于地球物理学的角动量守恒原理
  - 支持不同维度的气压厚度输入
  - 自动处理缺失值
  - 考虑地球半径和重力加速度
  - 使用余弦纬度加权
  - 适用于气候和大气环流分析

## 4. accumrun.f [📝 src](fortran/accumrun.f)
**功能**: 计算滑动窗口累积和(Running Accumulation)

**详细说明**:
- 包含一个子程序:
  - DACUMRUN: 计算指定窗口大小的滑动累积和
- 输入参数:
  - P: 输入时间序列数据数组(长度NP)
  - NP: 数据点总数
  - PMSG: 缺失值代码
  - NRUN: 滑动窗口大小
  - IOPT: 缺失值处理选项(0或1)
- 输出参数:
  - PAC: 累积和数组
  - NMSG: 遇到的缺失值总数
- 功能特点:
  - 前(NRUN-1)个值设为缺失值
  - 从第NRUN个点开始计算滑动累积和
  - 支持两种缺失值处理策略:
    - IOPT=0: 遇到缺失值时终止当前窗口计算，设结果为缺失值
    - IOPT=1: 遇到缺失值时计数但继续计算剩余有效值
  - 适用于降水累积、温度累积等时间序列分析

## 5. areaAve.f [📝 src](fortran/areaAve.f)
**功能**: 计算二维数据的加权区域平均值(使用分离权重)

**详细说明**:
- 包含一个子程序:
  - DWGTAREAAVE: 使用X和Y方向分离权重计算区域平均值
- 输入参数:
  - T: 2D数据数组(MX×NY)
  - WGTX: X方向权重数组(长度MX，如经度权重)  
  - WGTY: Y方向权重数组(长度NY，如纬度权重)
  - MX,NY: 数组维度(经度×纬度)
  - TMSG: 缺失值代码
  - IFLAG: 缺失值处理标志(0=忽略缺失值,1=遇到缺失值返回缺失值)
- 输出参数:
  - AVE: 加权区域平均值
- 功能特点:
  - 支持经纬度网格的面积加权
  - 权重为WGTX(i)×WGTY(j)的乘积
  - 灵活的缺失值处理策略
  - 适用于气候数据的空间平均计算

## 6. areaAve2.f [📝 src](fortran/areaAve2.f)  
**功能**: 计算二维数据的加权区域平均值(使用2D权重矩阵)

**详细说明**:
- 包含一个子程序:
  - DWGTAREAAVE2: 使用2D权重矩阵计算区域平均值
- 输入参数:
  - T: 2D数据数组(MX×NY)
  - WGT: 2D权重矩阵(MX×NY)
  - MX,NY: 数组维度  
  - TMSG: 缺失值代码
  - IFLAG: 缺失值处理标志
- 输出参数:
  - AVE: 加权区域平均值
- 功能特点:
  - 支持任意复杂的权重分布
  - 每个网格点可有独立权重值
  - 适用于不规则网格或特殊加权需求
  - 灵活的缺失值处理

## 7. areaRmse.f [📝 src](fortran/areaRmse.f)
**功能**: 计算两个二维数据的加权均方根误差(使用分离权重)

**详细说明**:
- 包含一个子程序:
  - DWGTAREARMSE: 使用X和Y方向分离权重计算RMSE
- 输入参数:
  - T,Q: 两个2D数据数组(MX×NY)，用于比较
  - WGTX: X方向权重数组(长度MX)
  - WGTY: Y方向权重数组(长度NY)  
  - MX,NY: 数组维度
  - TMSG,QMSG: T和Q的缺失值代码
  - IFLAG: 缺失值处理标志
- 输出参数:
  - RMSE: 加权均方根误差
- 功能特点:
  - 计算公式: $\text{RMSE} = \sqrt{\frac{\sum \text{wgt}(T-Q)^2}{\sum \text{wgt}}}$
  - 支持两个数据集的误差分析
  - 权重为$\text{WGTX}(i) \times \text{WGTY}(j)$的乘积
  - 适用于模型验证和数据对比

## 8. areaRmse2.f [📝 src](fortran/areaRmse2.f)  
**功能**: 计算两个二维数据的加权均方根误差(使用2D权重矩阵)

**详细说明**:
- 包含一个子程序:
  - DWGTAREARMSE2: 使用2D权重矩阵计算RMSE
- 输入参数:
  - T,Q: 两个2D数据数组(MX×NY)
  - WGT: 2D权重矩阵(MX×NY)
  - MX,NY: 数组维度
  - TMSG,QMSG: 缺失值代码  
  - IFLAG: 缺失值处理标志
- 输出参数:
  - RMSE: 加权均方根误差
- 功能特点:
  - 支持任意复杂的权重分布
  - 每个网格点可有独立权重
  - 适用于不规则网格的误差分析
  - 灵活处理缺失值

## 9. areaSum2.f [📝 src](fortran/areaSum2.f)
**功能**: 计算二维数据的加权区域总和(使用2D权重矩阵)

**详细说明**:
- 包含一个子程序:
  - DWGTAREASUM2: 使用2D权重矩阵计算区域加权总和
- 输入参数:
  - T: 2D数据数组(MX×NY)
  - WGT: 2D权重矩阵(MX×NY)
  - MX,NY: 数组维度
  - TMSG: 缺失值代码
  - IFLAG: 缺失值处理标志
- 输出参数:
  - SUM: 加权区域总和
- 功能特点:
  - 计算公式: $\text{SUM} = \sum_{i=1}^{MX} \sum_{j=1}^{NY} T(i,j) \times \text{WGT}(i,j)$
  - 支持任意复杂的权重分布
  - 每个网格点可有独立权重值
  - 灵活的缺失值处理策略:
    - IFLAG=0: 忽略缺失值，继续计算有效值的总和
    - IFLAG=1: 遇到缺失值时立即返回缺失值
  - 适用于不规则网格的面积积分计算
  - 灵活处理缺失值情况

## 10. area_conremap.f [📝 src](fortran/area_conremap.f)
**功能**: 保守插值的格点重映射(Conservative Remapping)

**详细说明**:
- 包含四个子程序:
  - CREMAPBIN: 主重映射程序，网格盒装箱方法
  - BINNING_GET_GLOBAL_LATS_WGTS: 获取全球高斯纬度和权重
  - BINNING_MAP_LATS: 将区域纬度映射到全球纬度
  - BINNING_MAP_EDGES: 计算网格盒边界
- 输入参数:
  - XX: 输入3D数据场(经度×纬度×层数)
  - CLAT,CLON: 输入网格的纬度和经度
  - CLATO,CLONO: 输出网格的纬度和经度  
  - BIN_FACTOR: 装箱面积扩展/收缩因子
  - NLAT,NLATO: 全球高斯纬度数
- 输出参数:
  - YY: 水平插值后的输出场
- 功能特点:
  - 支持规则和高斯网格
  - 保守插值算法，保持质量守恒
  - 支持垂直多层数据
  - 自动处理缺失值
  - 适用于气候模型数据重网格化

## 11. bfilter1.f [📝 src](fortran/bfilter1.f)  
**功能**: Butterworth带通滤波器

**详细说明**:
- 包含两个子程序:
  - BUTTFILT: 滤波器入口程序
  - BFILTER: 核心滤波器实现
- 输入参数:
  - XR: 输入时间序列(长度N)
  - M: Butterworth滤波器阶数
  - FCA,FCB: 低频和高频截止频率
  - DT: 时间间隔
  - MZER: 是否移除均值的标志
- 输出参数:
  - YR: 滤波后的时间序列
  - ER: 包络函数
  - IER: 错误码
- 功能特点:
  - 零相位Butterworth带通滤波
  - 级联一阶滤波器实现
  - 前向和反向滤波确保零相位
  - 复数时间序列处理
  - 稳定的滤波算法
  - 适用于地震学和信号处理

## 12. binavg.f [📝 src](fortran/binavg.f)
**功能**: 不规则分布数据的网格装箱平均

**详细说明**:
- 包含一个子程序:
  - BINDATAAVG: 将不规则地理数据进行网格化平均
- 输入参数:
  - ZLON,ZLAT: 数据点的经纬度坐标(NZ个点)
  - Z: 数据值数组(长度NZ)
  - GLON,GLAT: 目标网格的经纬度(MLON×NLAT)
  - ZMSG: 缺失值标识
  - IOPT: 选项参数
- 输出参数:
  - GBINKNT: 网格装箱结果(第1层为平均值，第2层为计数)
  - IER: 错误代码
- 功能特点:
  - 验证网格等间距性
  - 将数据点分配到对应网格单元
  - 计算每个网格单元内的平均值
  - 适用于地理数据重新格点化

## 13. bindata.f [📝 src](fortran/bindata.f)
**功能**: 不规则分布数据的网格装箱求和

**详细说明**:
- 包含一个子程序:
  - BINDATASUM3: 将不规则地理数据进行网格化求和
- 输入参数:
  - ZLON,ZLAT: 数据点坐标
  - Z: 数据值数组
  - GLON,GLAT: 目标网格坐标
  - ZMSG: 缺失值标识
- 输出参数:
  - GBIN: 网格装箱求和结果
  - GKNT: 每个网格单元的数据点计数
- 功能特点:
  - 执行地理数据网格装箱求和
  - 提供数据重新分布统计
  - 用于统计分析的预处理
  - 不计算平均值，只进行求和

## 14. calc_uh.f90 [📝 src](fortran/calc_uh.f90)
**功能**: 计算上升气流螺旋度(Updraft Helicity)

**详细说明**:
- 包含一个子程序:
  - DCALCUH: 计算UH诊断量
- 输入参数:
  - US,VS: u和v风速分量(nx×ny×nz)
  - W: 垂直速度(nx×ny×nzp1)
  - ZP: 高度数组
  - MAPFCT: 地图投影因子
  - UHMNHGT,UHMXHGT: 计算高度范围
- 输出参数:
  - UH: 上升气流螺旋度场(nx×ny)
- 功能特点:
  - 基于Kain等人2008年公式
  - 在指定高度范围内积分垂直涡度×垂直速度
  - 用于识别超级单体风暴和龙卷风
  - 使用OpenMP并行计算
  - 重要的气象诊断量

## 15. calcorc_dp.f [📝 src](fortran/calcorc_dp.f)
**功能**: 计算皮尔逊相关系数和协方差

**详细说明**:
- 包含两个子程序:
  - DCALCORC: 计算相关系数
  - DCALCOVC: 计算协方差
- 输入参数:
  - X,Y: 两个数据序列(长度NXY)
  - XAVE,XSTD: X序列的均值和标准差
  - XMSG,YMSG: 缺失值标识
- 输出参数:
  - CORC: 皮尔逊相关系数
  - COVC: 协方差
  - IER: 错误代码
- 功能特点:
  - 能够处理缺失值
  - 自动计算统计量
  - 标准化和未标准化计算
  - 统计分析基础工具

## 16. cancor.f [📝 src](fortran/cancor.f)
**功能**: 典型相关分析(Canonical Correlation Analysis)

**详细说明**:
- 包含多个子程序:
  - DCANCORXY: 主分析程序
  - DCANORX,DMINVX,DNROOTX,DEIGENX: 矩阵运算子程序
- 输入参数:
  - X: 独立变量矩阵(nobs×mx)
  - Y: 依赖变量矩阵(nobs×my)
  - NOBS: 观测数量
- 输出参数:
  - CANR: 典型相关系数
  - EVAL: 特征值
  - CHISQ: 卡方统计量
  - COEFXR,COEFYL: 典型系数
- 功能特点:
  - 识别两组变量间的线性关系模式
  - 包含完整矩阵运算套件
  - 多元统计分析高级技术
  - 用于降维和模式识别

## 17. cdf_dp.f [📝 src](fortran/cdf_dp.f)
**功能**: 多种概率分布的累积分布函数计算

**详细说明**:
- 包含多个子程序:
  - DCDFBINP,DCDFBINS,DCDFBINXN: 二项分布CDF
  - DCDFGAMP,DCDFGAMX: 伽马分布CDF
  - DCDFNORP,DCDFNORX: 正态分布CDF
  - DCDFCHIP: 卡方分布CDF
- 输入参数: 各种分布的参数(概率、自由度、形状参数等)
- 输出参数: 对应的CDF值或逆函数值
- 功能特点:
  - 支持多种概率分布
  - 提供正向和逆向计算
  - 完善的错误处理
  - 概率统计和假设检验核心工具

## 18. cdft_dp.f [📝 src](fortran/cdft_dp.f)  
**功能**: t分布的累积分布函数计算

**详细说明**:
- 包含三个子程序:
  - DCDFCDFTT,DCDFCDFTP: t分布CDF计算
  - DCDFTERR: 错误处理
- 输入参数:
  - P: 概率值数组
  - T: t统计量数组
  - DF: 自由度数组
- 输出参数:
  - 根据计算方向输出t值或概率值
  - IER: 错误代码
- 功能特点:
  - 专门处理t分布计算
  - 提供双向计算功能
  - 完善的错误处理机制
  - t检验和置信区间计算基础工具

## 19. cfft_driver.f [📝 src](fortran/cfft_driver.f)
**功能**: 复数FFT驱动程序

**详细说明**:
- 包含四个子程序:
  - CFFTFDRIVER: 复数前向FFT驱动
  - CFFTBDRIVER: 复数逆向FFT驱动
  - FRQCFFT: 生成FFT频率
  - CFFTFFRQREORDER: 重排序频率
- 输入参数: 复数实部虚部数组、工作数组
- 输出参数: FFT系数、频率数组
- 功能特点:
  - 为NCL提供复数FFT高级接口
  - 包括前向变换、逆向变换功能
  - 支持频率生成和重排序
  - 适用于时间序列频谱分析和信号处理

## 20. cfftb.f [📝 src](fortran/cfftb.f)
**功能**: 复数FFT逆变换

**详细说明**:
- 包含一个子程序:
  - CFFTB: 复数离散傅里叶逆变换
- 输入参数: 复数序列C、工作数组WSAVE
- 输出参数: 逆变换后的复数序列
- 功能特点:
  - 执行傅里叶合成
  - 从傅里叶系数重构复数周期序列
  - 频域到时域转换
  - 用于信号重构和滤波后恢复

## 21. cfftb1.f [📝 src](fortran/cfftb1.f)
**功能**: 复数FFT逆变换核心算法

**详细说明**:
- 包含一个子程序:
  - CFFTB1: 复数FFT逆变换核心实现
- 输入参数: 序列长度N、复数数组C、工作数组
- 输出参数: 变换结果
- 功能特点:
  - 底层逆变换实现
  - 处理不同因子分解的FFT计算
  - 支持2、3、4、5等基本因子
  - FFT算法内部核心计算

## 22. cfftf.f [📝 src](fortran/cfftf.f)
**功能**: 复数FFT前向变换

**详细说明**:
- 包含一个子程序:
  - CFFTF: 复数离散傅里叶前向变换
- 输入参数: 复数序列C、工作数组WSAVE
- 输出参数: 傅里叶系数
- 功能特点:
  - 执行傅里叶分析
  - 计算复数周期序列的傅里叶系数
  - 时域到频域转换
  - 用于频谱分析和信号分解

## 23. cfftf1.f [📝 src](fortran/cfftf1.f)
**功能**: 复数FFT前向变换核心算法

**详细说明**:
- 包含一个子程序:
  - CFFTF1: 复数FFT前向变换核心实现
- 输入参数: 序列长度N、复数数组C、工作数组
- 输出参数: 变换结果
- 功能特点:
  - 底层前向变换实现
  - 采用混合基算法
  - 处理不同长度的序列
  - 高效前向变换处理

## 24. cffti.f [📝 src](fortran/cffti.f)
**功能**: 复数FFT初始化

**详细说明**:
- 包含一个子程序:
  - CFFTI: 复数FFT初始化程序
- 输入参数: 序列长度N
- 输出参数: 工作数组WSAVE
- 功能特点:
  - 初始化FFT计算所需工作数组
  - 计算并存储素因子分解
  - 生成三角函数表
  - FFT预处理和工作数组准备

## 25. cffti1.f [📝 src](fortran/cffti1.f)
**功能**: 复数FFT初始化核心算法

**详细说明**:
- 包含一个子程序:
  - CFFTI1: 复数FFT初始化核心实现
- 输入参数: 序列长度N
- 输出参数: 工作数组WA、因子数组IFAC
- 功能特点:
  - 执行FFT实际初始化工作
  - 素因子分解
  - 三角函数计算和存储
  - 优化重复变换的预计算

## 26. chiin_dp.f [📝 src](fortran/chiin_dp.f)
**功能**: 卡方分布逆函数计算

**详细说明**:
- 包含一个子程序:
  - CHISUB: 卡方分布逆函数接口
- 输入参数: 累积概率P、自由度DF
- 输出参数: 卡方值CHI
- 功能特点:
  - 计算卡方分布的逆函数
  - 给定概率和自由度求卡方统计量
  - 适用于统计假设检验
  - 置信区间计算和概率分析

## 27. cntrFiniteDiff.f [📝 src](fortran/cntrFiniteDiff.f)
**功能**: 中心有限差分计算

**详细说明**:
- 包含一个子程序:
  - DCFINDIF: 中心有限差分计算
- 输入参数: 变量数组Q、坐标数组R、边界条件参数
- 输出参数: 导数数组DQDR
- 功能特点:
  - 使用中心有限差分方法
  - 计算数值导数
  - 支持循环和非循环边界条件
  - 用于数值微分和梯度计算

## 28. conservLinInt_2d_ncl.f [📝 src](fortran/conservLinInt_2d_ncl.f)
**功能**: 二维保守线性插值

**详细说明**:
- 包含一个子程序:
  - AREALININT2DA: 二维面积加权插值
- 输入参数: 输入网格坐标、数据、权重
- 输出参数: 插值后的网格数据
- 功能特点:
  - 实现二维保守线性插值
  - 面积加权方法
  - 保持数据的积分守恒
  - 适用于网格重映射和气候数据处理

## 29. covcorm_driver.f [📝 src](fortran/covcorm_driver.f)
**功能**: 协方差相关矩阵计算驱动

**详细说明**:
- 包含两个子程序:
  - DCOVCORMSSM: 对称存储模式协方差矩阵
  - DCOVCORM: 二维协方差矩阵
- 输入参数: 数据矩阵、缺失值标识、选项
- 输出参数: 协方差或相关矩阵
- 功能特点:
  - 计算协方差矩阵或相关矩阵
  - 支持对称存储和全矩阵存储
  - 多变量统计分析工具
  - 用于主成分分析和相关性分析

## 30. covcorm_xy_matrix_ncl.f [📝 src](fortran/covcorm_xy_matrix_ncl.f)
**功能**: 交叉协方差矩阵计算

**详细说明**:
- 包含一个子程序:
  - DCOVARXY: 双变量协方差计算
- 输入参数: 变量数组X、Y、滞后参数
- 输出参数: 交叉协方差/相关矩阵
- 功能特点:
  - 计算两个多维变量间的交叉协方差
  - 支持时间滞后分析
  - 交叉相关分析工具
  - 适用于滞后相关研究和时间序列分析

## 31. ctwrap.f [📝 src](fortran/ctwrap.f)
**功能**: 三维网格可视化系统

**详细说明**:
- 包含多个子程序: CTDRIVER, DFCLRS, DRWMCL, DCFOCB, MKMESH, DSJE02等
- 输入参数: 三角网格数据、绘图参数、工作数组
- 输出参数: 各种可视化图形和三角网格
- 功能特点:
  - 复杂的三维网格可视化系统
  - 从SEAM格网生成三角网格
  - 支持多种可视化方式
  - 用于地球表面网格数据的三维可视化、等高线绘制、填色图制作

## 32. cz2ccm_dp.f [📝 src](fortran/cz2ccm_dp.f)
**功能**: CCM2混合坐标系统地球物理高度计算

**详细说明**:
- 包含两个子程序: DCZ2CCM, DCZ2
- 输入参数: 地表压力PS、地表位势PHIS、虚温TV、混合坐标系数
- 输出参数: 中间层位势高度Z2
- 功能特点:
  - 基于CCM2混合坐标系统
  - 计算地球物理高度场
  - 用于大气模式中的垂直坐标转换和气象数据处理

## 33. d1mach_alpha.f [📝 src](fortran/d1mach_alpha.f)
**功能**: Alpha平台双精度浮点常数

**详细说明**:
- 包含一个子程序: D1MACH
- 输入参数: 整数参数I (1-5)
- 输出参数: Alpha架构机器相关的双精度浮点常数
- 功能特点:
  - DEC Alpha处理器优化的浮点常数
  - IEEE 754标准兼容
  - 高精度数值计算支持

## 34. d1mach_linux.f [📝 src](fortran/d1mach_linux.f)
**功能**: Linux平台双精度浮点常数

**详细说明**:
- 包含一个子程序: D1MACH
- 输入参数: 整数参数I (1-5)
- 输出参数: Linux系统机器相关的双精度浮点常数
- 功能特点:
  - x86/x64架构优化
  - GNU编译器兼容
  - 跨Linux发行版支持

## 35. d1mach_sun.f [📝 src](fortran/d1mach_sun.f)
**功能**: Sun Solaris平台双精度浮点常数

**详细说明**:
- 包含一个子程序: D1MACH
- 输入参数: 整数参数I (1-5)
- 输出参数: Sun SPARC架构机器相关的双精度浮点常数
- 功能特点:
  - SPARC处理器优化
  - Sun编译器兼容
  - Solaris系统集成

## 36. dbetai.f [📝 src](fortran/dbetai.f)
**功能**: 不完全Beta函数计算

**详细说明**:
- 包含主程序: BETAINC, DBETAISLATEC及多个SLATEC辅助函数
- 输入参数: 积分上限x、beta分布参数a,b
- 输出参数: 不完全Beta函数值result
- 功能特点:
  - 计算不完全Beta函数
  - 包含完整的SLATEC数学库支持
  - 用于统计分析、概率计算、假设检验

## 37. dcolconv.f [📝 src](fortran/dcolconv.f)
**功能**: 颜色空间转换库

**详细说明**:
- 包含多个子程序: DHLSRGB, DHSVRGB, DRGBHLS, DRGBHSV, DRGBYIQ, DYIQRGB等
- 输入参数: 各种颜色空间的颜色值
- 输出参数: 转换后的颜色值
- 功能特点:
  - 完整的颜色空间转换库
  - 支持RGB、HLS、HSV、YIQ之间互相转换
  - 用于图形可视化、颜色处理、图像处理

## 38. det_code42.f [📝 src](fortran/det_code42.f)
**功能**: 矩阵行列式计算

**详细说明**:
- 包含三个子程序: DETCALC, DTRM, ELGS
- 输入参数: n×n矩阵a
- 输出参数: 矩阵行列式值det
- 功能特点:
  - 使用部分选主元高斯消去法
  - 计算矩阵行列式
  - 用于线性代数计算、数值分析、求解线性方程组

## 39. dewpt_trh.f [📝 src](fortran/dewpt_trh.f)
**功能**: 露点温度计算

**详细说明**:
- 包含一个子程序: DDEWTEMP
- 输入参数: 温度TK[K]、相对湿度RH[%]
- 输出参数: 露点温度TDK[K]
- 功能特点:
  - 根据温度和相对湿度计算露点温度
  - 用于气象分析、湿度计算、大气物理研究

## 40. dgeevx_interface.f [📝 src](fortran/dgeevx_interface.f)
**功能**: LAPACK特征值计算接口

**详细说明**:
- 包含一个子程序: DGEEVXINT
- 输入参数: 矩阵A、控制参数、工作数组
- 输出参数: 特征值WR,WI和特征向量EVLR
- 功能特点:
  - LAPACK DGEEVX函数的接口
  - 计算矩阵特征值和特征向量
  - 用于数值线性代数、主成分分析、动力系统分析

## 41. dimavgwgt.f [📝 src](fortran/dimavgwgt.f)
**功能**: 一维加权平均值计算

**详细说明**:
- 包含一个子程序: DIMAVGWGT
- 输入参数: 数组x、权重wgt、控制参数iopt
- 输出参数: 加权平均值xavg
- 功能特点:
  - 计算一维数组的加权平均值
  - 支持缺失值处理
  - 用于统计分析、数据平均化、带权重的数据处理

## 42. dimsumwgt.f [📝 src](fortran/dimsumwgt.f)
**功能**: 一维加权和计算

**详细说明**:
- 包含一个子程序: DIMSUMWGT
- 输入参数: 数组x、权重wgt、控制参数iopt
- 输出参数: 加权和sumx
- 功能特点:
  - 计算一维数组的加权和
  - 支持缺失值处理
  - 用于统计分析、数据求和、带权重的数据聚合

## 43. distribfit.f [📝 src](fortran/distribfit.f)
**功能**: Weibull分布参数拟合

**详细说明**:
- 包含两个子程序: WEIBFIT, QNORMALQ
- 输入参数: 数据数组xin、最小样本数nmin、置信水平level
- 输出参数: Weibull分布参数weib(6)
- 功能特点:
  - 使用最大似然法拟合双参数Weibull分布
  - 用于可靠性分析、生存分析、极值统计、故障分析

## 44. dmapgci.f [📝 src](fortran/dmapgci.f)
**功能**: 地球表面大圆插值

**详细说明**:
- 包含两个子程序: DMAPGCI, DGCDIST
- 输入参数: 两点的经纬度、插值点数量
- 输出参数: 大圆路径上的插值点坐标
- 功能特点:
  - 在地球表面两点间的大圆路径上进行插值
  - 计算大圆距离
  - 用于地理信息系统、航行路线规划、地理数据插值

## 45. dpsort.f [📝 src](fortran/dpsort.f)
**功能**: 双精度被动排序

**详细说明**:
- 包含两个子程序: DPSORTDRIVER, DPSORT
- 输入参数: 数组dx、排序控制参数kflag
- 输出参数: 排序索引数组iperm
- 功能特点:
  - 被动排序算法，使用修改的快速排序
  - 产生排列向量而不改变原数组
  - 用于数据排序、索引生成、数据重排序

## 46. dpsort_large.f [📝 src](fortran/dpsort_large.f)
**功能**: 大数组被动排序

**详细说明**:
- 包含两个子程序: DPSORTLARGEDRIVER, DPSORTLARGE
- 输入参数: 大数组dx、排序控制参数kflag
- 输出参数: 排序索引数组iperm (64位整数)
- 功能特点:
  - 大数组的被动排序算法
  - 支持64位索引
  - 用于大数据排序、超大数组处理、高性能计算

## 47. dtrend_dp.f [📝 src](fortran/dtrend_dp.f)
**功能**: 时间序列趋势移除

**详细说明**:
- 包含两个子程序: DDTRNDX, DDTRNDMSG
- 输入参数: 时间序列数据、去趋势选项
- 输出参数: 去趋势后的数据、趋势系数
- 功能特点:
  - 时间序列的趋势移除
  - 支持线性和二次趋势
  - 允许缺失值，用于时间序列分析、气候数据处理、趋势分析

## 48. dtrend_lsq_msg_dp.f [📝 src](fortran/dtrend_lsq_msg_dp.f)
**功能**: 最小二乘二次趋势移除

**详细说明**:
- 包含一个子程序: QDTRNDMSG
- 输入参数: 时间序列Y、缺失值标志、控制选项
- 输出参数: 去除二次趋势后的序列、趋势系数
- 功能特点:
  - 使用最小二乘法进行二次趋势移除
  - 支持缺失值
  - 用于时间序列预处理、非线性趋势分析、数据去噪

## 49. dz_height.f [📝 src](fortran/dz_height.f)
**功能**: 大气层厚度计算

**详细说明**:
- 包含三个子程序: DZHGTDRV, DZHEIGHT, DZMONO
- 输入参数: 高度场z、地表高度zsfc、顶部高度ztop
- 输出参数: 各层厚度dz
- 功能特点:
  - 计算大气层厚度
  - 基于给定的高度场和边界条件
  - 用于大气模式、垂直坐标处理、层厚计算

## 50. eof2data.f [📝 src](fortran/eof2data.f)
**功能**: EOF分析数据重构

**详细说明**:
- 包含一个子程序:
  - DEOF2DATA: 从EOF模态和时间系数重构数据场
- 输入参数:
  - NEVAL,NPTS,NTIM: 模态数、空间点数、时间点数
  - EV: EOF特征向量(NPTS×NEVAL)
  - EVTS: 主成分时间系数(NTIM×NEVAL)
  - XMSG: 缺失值代码
- 输出参数:
  - X: 重构的时空数据场(NTIM×NPTS)
- 功能特点:
  - 重构公式: $X(t,s) = \sum_k[\text{EVTS}(t,k) \times \text{EV}(s,k)]$
  - 支持部分模态重构(使用前NEVAL个模态)
  - 自动处理缺失值情况
  - 适用于EOF降维后的数据恢复和重构验证

## 51. eof_scripps.f [📝 src](fortran/eof_scripps.f)
**功能**: 经验正交函数(EOF)分析

**详细说明**:
- 包含多个子程序:
  - DEOF11: 计算一个数据集的EOF分析，使用协方差矩阵
  - DEOF22: 计算两个数据集间的交叉EOF分析
  - DEOF: 主要EOF计算程序
  - DEOFCOVCOR: 协方差/相关矩阵计算支持
  - DSYMMEIG: 对称矩阵特征值分解
- 输入参数:
  - X: 时空数据场(time×space)
  - NTIME,NSPACE: 时间和空间维度
  - NEVAL: 需要计算的特征向量个数
  - XMSG: 缺失值代码
  - IOPT: 计算选项(协方差=0，相关矩阵=1)
- 输出参数:
  - EVAL: 特征值(方差贡献)
  - EVEC: 特征向量(EOF模态)
  - PCF: 主成分时间系数
  - PCVAR: 各模态方差贡献百分比
- 功能特点:
  - 基于经典的EOF分解理论
  - 支持缺失值的处理
  - 提供协方差和相关矩阵两种分析模式
  - 计算方差贡献和主成分
  - 适用于气候模态分析、降维、模式识别

## 52. eqthecalc.f90 [📝 src](fortran/eqthecalc.f90)
**功能**: 等效位温(Equivalent Potential Temperature)计算

**详细说明**:
- 包含一个主程序:
  - DEQTHECALC: WRF模式等效位温计算程序
- 输入参数:
  - QVP: 水汽混合比(g/kg)，3D数组(miy×mjx×mkzh)
  - TMK: 温度(K)，3D数组(miy×mjx×mkzh)
  - PRS: 全压力(hPa)，3D数组(miy×mjx×mkzh)
  - MIY,MJX,MKZH: 数组维度
- 输出参数:
  - ETH: 等效位温(K)，3D数组(miy×mjx×mkzh)
- 功能特点:
  - 使用WRF物理常数模块(EPS,GAMMA等)
  - 计算抬升凝结高度TLCL
  - 等效位温公式: $\theta_e = T \times \left(\frac{1000}{p}\right)^{\gamma} \times \exp[\text{热力学项}]$
  - 使用OpenMP并行化(COLLAPSE(3))
  - 适用于大气对流分析、不稳定度诊断

## 53. erf_sub.f [📝 src](fortran/erf_sub.f)
**功能**: 误差函数计算简单接口

**详细说明**:
- 包含两个子程序:
  - DERRF: 调用标准erf函数的包装器
  - DERRCF: 调用erfc1函数的包装器
- 输入参数:
  - DERRF: X(标量自变量)
  - DERRCF: IND(整数参数), X(自变量)  
- 输出参数:
  - RESULT: 误差函数计算结果
- 功能特点:
  - 提供NCL对标准误差函数的简单接口
  - $\text{erf}(x) = \frac{2}{\sqrt{\pi}}\int_0^x e^{-t^2}dt$  
  - $\text{erfc}(x) = 1 - \text{erf}(x)$
  - 单点计算而非向量化
  - 用于概率计算、统计分析

## 54. evssm_lapack.f [📝 src](fortran/evssm_lapack.f)
**功能**: 对称矩阵特征值/特征向量计算

**详细说明**:
- 包含一个子程序:
  - EVSSM: 使用LAPACK DSPEVX的特征值求解器
- 输入参数:
  - NVAR: 矩阵阶数
  - LSSM: 对称存储矩阵大小
  - CSSM: 对称矩阵(紧致存储)
  - NEVAL: 需要计算的特征值个数
  - IOPT: 排序选项(0=降序,1=升序)
- 输出参数:
  - WEVAL: 特征值数组
  - WEVEC: 特征向量矩阵(NVAR×NEVAL)
  - INFO: LAPACK返回状态
- 功能特点:
  - 调用LAPACK DSPEVX例程
  - 计算最大的NEVAL个特征值和对应特征向量
  - 支持升序/降序排列选择
  - 适用于主成分分析、协方差矩阵分解

## 55. exptapersh.f [📝 src](fortran/exptapersh.f)
**功能**: 球谐系数指数锥形滤波

**详细说明**:
- 包含两个子程序:
  - DEXPTAPER: 计算指数锥形滤波权重
  - DEXPTAPERSH: 对球谐系数应用滤波
- 输入参数:
  - DEXPTAPER: N0(截止参数), R(指数), NLAT(纬度数)
  - DEXPTAPERSH: AB,BB(球谐系数), N0, R
- 输出参数:
  - DEXPTAPER: S(滤波权重数组)
  - DEXPTAPERSH: 滤波后的AB,BB系数
- 功能特点:
  - 滤波权重: $S(n) = \exp\left(-\left(\frac{n(n+1)}{N_0(N_0+1)}\right)^R\right)$
  - 对球谐系数AB(m,n)和BB(m,n)应用权重S(n)
  - 抑制高阶球谐波，保持低阶分量
  - 用于全球数据的频谱平滑和去噪

## 56. ezfft_dp.f [📝 src](fortran/ezfft_dp.f)
**功能**: FFTPACK双精度FFT子集

**详细说明**:
- 包含多个子程序:
  - DEZFFTF: 简化实数前向FFT(输出Fourier系数)
  - DEZFFTB: 简化实数逆向FFT(从Fourier系数重构)
  - DEZFFTI: EZ-FFT初始化程序
  - DRFFTF/DRFFTB: 标准实数FFT前向/逆向变换
  - DRFFTI: 实数FFT初始化
- 输入参数:
  - N: 序列长度
  - R: 输入实数序列
  - WSAVE: 预计算的工作数组
- 输出参数:
  - AZERO: 直流分量
  - A,B: 余弦和正弦系数数组
- 功能特点:
  - FFTPACK库的双精度版本
  - 支持混合基FFT算法(2,3,4,5因子分解)
  - EZFFT提供用户友好的系数格式
  - 高效的实数FFT专用算法
  - 适用于实数时间序列的频谱分析

## 57. f2foshv.f [📝 src](fortran/f2foshv.f)
**功能**: 矢量场固定网格到固定偏移网格变换

**详细说明**:
- 包含一个主程序:
  - DF2FOSHV: 调用SPHEREPACK的VSHIFTI和VSHIFTE例程
- 输入参数:
  - UOFF,VOFF: 固定偏移网格上的u,v分量(ILON×JLAT)
  - UREG,VREG: 固定网格上的u,v分量(ILON×JLAT1)
  - ILON,JLAT: 网格维度
  - WORK,WSAVE: 工作数组
- 输出参数:
  - IER: 错误状态码
- 功能特点:
  - 基于SPHEREPACK库的球面矢量场变换
  - 处理球面坐标系的矢量分量正确变换
  - 保持矢量场的物理一致性
  - 用于大气海洋模式网格变换

## 58. fill_msg_grid.f [📝 src](fortran/fill_msg_grid.f)
**功能**: 泊松方程法填充网格缺失值

**详细说明**:
- 包含两个子程序:
  - POISXY1: 缺失值填充的入口程序
  - POISXY2: 泊松方程松弛迭代求解器
- 输入参数:
  - XIO: 含缺失值的2D数组(MX×NY)
  - XMSG: 缺失值标识
  - GTYPE: 边界条件(0=非循环, 1=x方向循环)
  - NSCAN: 最大迭代次数
  - EPSX,RELC: 收敛标准和松弛常数
  - GUESS: 初始猜测(0=零值, 1=纬向平均)
- 输出参数:
  - 填充后的XIO数组, MSCAN: 实际迭代次数
- 功能特点:
  - 使用拉普拉斯方程: $\nabla^2 Z = 0$
  - 四点平均松弛迭代法
  - 支持循环和非循环边界条件
  - 适用于空间数据插值和缺失值填补

## 59. filtrx.f [📝 src](fortran/filtrx.f)
**功能**: Lanczos和高斯数字滤波器设计

**详细说明**:
- 包含三个子程序:
  - DFILTRQ: Lanczos滤波器权重和响应函数计算
  - DFILWTQ: 滤波器权重计算和归一化
  - FILWGTNORMAL: 高斯(正态)滤波器权重生成
- 输入参数:
  - DFILTRQ: NWT(权重数), FCA/FCB(截止频率), IHP(滤波类型)
  - FILWGTNORMAL: NWT, SIGMA(标准差), ITYPE
- 输出参数:
  - WT: 归一化滤波器权重
  - FREQ,RESP: 频率和响应数组
- 功能特点:
  - 支持低通(IHP=0)、高通(IHP=1)、带通(IHP=2)滤波
  - Lanczos窗减少Gibbs现象
  - 高斯滤波器基于正态分布
  - 提供完整的频率响应分析
  - 适用于时间序列滤波和信号处理

## 60. finfo.f [📝 src](fortran/finfo.f)
**功能**: 傅里叶级数信息分析

**详细说明**:
- 包含两个子程序:
  - DFOURINFO: NCL接口的傅里叶信息分析
  - DFINFO: 核心傅里叶级数分析程序
- 输入参数:
  - X: 输入时间序列(NPTS个点)
  - NPTS: 数据点数
  - NHRET: 需要返回的谐波数
  - SCLPHA: 相位比例因子
- 输出参数:
  - AMP: 各谐波振幅
  - PHA: 各谐波相位(经SCLPHA缩放)
  - PCV: 各谐波方差贡献百分比
- 功能特点:
  - 调用DEZFFTF进行傅里叶分解
  - 振幅: AMP(n) = $\sqrt{A(n)^2 + B(n)^2}$
  - 相位: PHASE(n) = atan2(B(n),A(n))×调整因子
  - 方差贡献: PCV(n) = $(0.5 \times \text{AMP}(n)^2/\text{总方差}) \times 100\%$
  - 适用于周期信号分析、谐波检测

## 61. fluxEddyTav_dp.f [📝 src](fortran/fluxEddyTav_dp.f)
**功能**: 涡动通量计算(双精度函数)

**详细说明**:
- 包含一个函数:
  - DFLXEDY: 计算两个时间序列的涡动通量
- 输入参数:
  - X,Y: 两个时间序列数组(NPTS)
  - NPTS: 数据点数
  - XMSG: 缺失值标识
- 输出参数:
  - DFLXEDY: 返回涡动通量值x'y'
  - IER: 错误状态码
- 功能特点:
  - 雷诺分解: $x = \bar{X} + x'$, $y = \bar{Y} + y'$
  - 涡动通量: x'y' = x̄ȳ - X̄Ȳ
  - 内部计算相关系数: $CC = \frac{x'y'}{\sqrt{x'x' \times y'y'}}$
  - 自动处理缺失值
  - 适用于大气湍流分析、热通量、动量通量计算

## 62. fo2fs_dp.f [📝 src](fortran/fo2fs_dp.f)
**功能**: 标量场固定偏移网格到固定网格变换

**详细说明**:
- 包含一个主程序:
  - DFO2F: 调用SPHEREPACK的DSSHIFTI和DSSHIFTE例程
- 输入参数:
  - GOFF: 固定偏移网格标量场(ILON×JLAT)
  - GREG: 固定网格标量场(ILON×JLAT1)
  - ILON,JLAT: 网格维度
  - WORK,WSAVE: 工作数组
- 输出参数:
  - IER: 错误状态码
- 功能特点:
  - 基于SPHEREPACK库的球面标量场变换
  - 处理球面坐标系的标量场正确变换
  - 用于标量场在不同网格间的转换
  - 大气海洋模式数据后处理

## 63. gamfitd3.f [📝 src](fortran/gamfitd3.f)
**功能**: 不完全伽马分布参数估计

**详细说明**:
- 包含一个主程序:
  - GAMFITD3: 估计伽马分布参数的最大似然方法
- 输入参数:
  - DATARR: 数据数组(N)
  - N: 数据点数
  - DMSG: 缺失值标识
  - PCRIT: 有效数据百分比临界值
  - INVSCALE: 是否返回1/$\beta$(0=否,1=是)
- 输出参数:
  - ALPHA: $\alpha$参数(log-space参数)
  - SCALE: $\beta$参数(尺度参数)
  - SHAPE: $\eta$参数(形状参数)
  - PZERO: 零值概率
  - IER: 状态码
- 功能特点:
  - 使用Thom(1958)最大似然估计公式
  - $\alpha = \log(\bar{x}) - (\sum\log(x))/n$
  - $\eta = (1+\sqrt{1+4\alpha/3})/(4\alpha)$
  - 处理零值和缺失值
  - 适用于降水分布拟合、极值分析

## 64. gamma_interface.f [📝 src](fortran/gamma_interface.f)
**功能**: SLATEC伽马函数向量化接口

**详细说明**:
- 包含一个主程序:
  - GAMMACOMPLETE: 调用SLATEC DGAMMASLATEC函数
- 输入参数:
  - NX: 数组长度
  - XIN: 输入自变量数组
  - HAS_MSG: 是否有缺失值(0=无,1=有)
  - XMSG: 缺失值标识
- 输出参数:
  - XOUT: 伽马函数值数组
- 功能特点:
  - 批量计算伽马函数值
  - 支持缺失值处理
  - 调用高精度SLATEC数学库
  - $\Gamma(x) = \int_0^\infty t^{x-1}e^{-t}dt$
  - 适用于概率分布计算、统计分析

## 65. gaus.f [📝 src](fortran/gaus.f)
**功能**: SPHEREPACK高斯求积点计算

**详细说明**:
- 包含四个子程序:
  - GAQDNCL: 主控程序，计算高斯纬度和权重
  - GAQDNCL1: 设置三对角矩阵并求解特征值问题
  - DRSTNCL: 实对称三对角矩阵特征值计算(EISPACK修改版)
  - DINTQLNCL: QL算法特征值迭代求解
  - DPYTHANCL: 计算$\sqrt{a^2+b^2}$避免溢出
- 输入参数:
  - NLAT: 高斯纬度点数
  - WORK: 工作数组(长度$\geq$4×NLAT×(NLAT+1)+2)
  - LWORK: 工作数组长度
- 输出参数:
  - THETA: 高斯纬度点(弧度，0到$\pi$)
  - WTS: 对应的高斯权重
  - IERROR: 错误代码
- 功能特点:
  - 基于勒让德多项式递推关系构造三对角矩阵
  - 特征值即为高斯点，第一个特征向量分量平方×2即为权重
  - 使用QL算法求解对称三对角矩阵特征值问题
  - 结果按弧度递增排列
  - 用于球面调和分析的高精度数值积分

## 66. gausLobat.f [📝 src](fortran/gausLobat.f)
**功能**: 高斯-洛巴托配置点和权重计算

**详细说明**:
- 包含多个子程序:
  - GAUSLOBAT: NCL接口，返回[-90,90]度的高斯-洛巴托纬度
  - FINDGLW: 根据给定纬度计算对应的高斯-洛巴托权重
  - GAUSSLOBATTO: 核心计算程序，求解高斯-洛巴托点和权重
  - NEWTONRAPHSON1/2: Newton-Raphson和二分法混合迭代求根
  - JACOBF: 计算勒让德多项式及其一、二阶导数
- 输入参数:
  - NPTS: 配置点数
  - XGLAT: 输入纬度数组(用于FINDGLW)
- 输出参数:
  - XGL: 高斯-洛巴托配置点(度或[-1,1])
  - WEIGHT: 对应的权重
- 功能特点:
  - 基于勒让德多项式 $P_N(x)$ 和其导数 $P'_N(x)$ 的零点
  - 端点包含在配置中：$x_1 = -1$，$x_{N+1} = 1$
  - 权重公式：$w_i = \dfrac{2}{N(N+1)\left[P_N(x_i)\right]^2}$
  - 高精度 Newton-Raphson 迭代求解（容差 $1\times10^{-15}$）
  - 用于谱方法、高精度数值积分

## 67. gendat.f [📝 src](fortran/gendat.f)
**功能**: 二维测试数据生成器

**详细说明**:
- 包含两个子程序:
  - DGENDAT: 生成具有指定特征的2D测试数据场
  - DFRAN: 简单的随机数生成器
- 输入参数:
  - DATA: 输出数据数组(IDIM×N)
  - M,N: 实际使用的数组维度
  - MLOW,MHGH: 低值和高值中心的数量(1-25)
  - DLOW,DHGH: 数据的最小和最大值
  - ISEED: 随机数种子(0-100)
- 输出参数:
  - 填充的DATA数组
- 功能特点:
  - 使用指数函数和生成数据场: $\sum \exp(-\text{距离}^2)$
  - 随机分布的高值和低值中心
  - 精确控制最小值和最大值
  - 100个预设随机数序列
  - 适用于图形程序测试、算法验证

## 68. hydro_dp.f [📝 src](fortran/hydro_dp.f)
**功能**: 静力学高度计算

**详细说明**:
- 包含一个子程序:
  - DHYDRO: 使用静力学方程计算位势高度
- 输入参数:
  - P: 气压数组(Pa，NLVL层)
  - TKV: 虚温数组(K，NLVL层)
  - ZSFC: 地面位势高度(gpm)
  - NLVL: 垂直层数
- 输出参数:
  - ZH: 各层位势高度(gpm)
  - IER: 错误代码
- 功能特点:
  - 静力学方程: dZ = -(RT/g)d(ln p)
  - 积分形式: $Z(k) = Z(k-1) + \frac{R}{g} \times \bar{T} \times \ln\left(\frac{p(k-1)}{p(k)}\right)$
  - T̄使用对数加权平均: T̄ = (T₁ln p₁ + T₂ln p₂)/ln(p₁p₂)
  - 使用WMO标准重力加速度: $g = 9.80665$ m/s$^2$
  - 干空气气体常数: $R = 287.04$ J/(kg·K)
  - 适用于大气物理高度计算

## 69. hyi2hyo.f [📝 src](fortran/hyi2hyo.f)
**功能**: 混合坐标系垂直插值

**详细说明**:
- 包含一个子程序:
  - DHYI2HYOB: 在不同混合坐标层间插值变量
- 输入参数:
  - P0: 基准压力(Pa)
  - HYAI,HYBI: 输入混合坐标系数(KLEVI层)
  - HYAO,HYBO: 输出混合坐标系数(KLEVO层)
  - PSFC: 地面压力(MLON×NLAT)
  - XI: 输入变量(MLON×NLAT×KLEVI)
  - INTFLG: 外推标志(0=设缺失值,1=用最近值)
  - XMSG: 缺失值标识
- 输出参数:
  - XO: 输出变量(MLON×NLAT×KLEVO)
  - MSGFLG: 缺失值存在标志
- 功能特点:
  - 混合坐标压力: $p(k) = \text{hya}(k) \times p_0 + \text{hyb}(k) \times p_{\text{sfc}}$
  - 对数线性插值: ln p 空间内的线性插值
  - 处理超出范围的情况(外推或设缺失值)
  - 适用于不同大气模式垂直坐标转换

## 70. index77_dp.f [📝 src](fortran/index77_dp.f)
**功能**: 气候指数计算(如南方涛动指数)

**详细说明**:
- 包含两个子程序:
  - DINDX77: 简化接口，管理工作数组分配
  - DINDX77X: 核心计算程序
- 输入参数:
  - X,Y: 两个站点的月度数据(NMOS×NYRS)
  - NMOS,NYRS: 月数和年数(通常NMOS=12)
  - XMSG: 缺失值标识
  - IPRNT: 打印标志
- 输出参数:
  - ZI: 气候指数(NMOS×NYRS)
  - ZNI: 噪声指数(NMOS×NYRS)
- 功能特点:
  - 基于Trenberth(1984)方法计算南方涛动指数
  - 计算各月气候均值和标准差
  - 总体异常标准差标准化: $\text{ZI} = \frac{X'}{\sigma_X} - \frac{Y'}{\sigma_Y}$
  - 噪声指数: $\text{ZNI} = \frac{X'}{\sigma_X} + \frac{Y'}{\sigma_Y}$
  - 适用于ENSO、AO等气候指数计算

## 71. int2p_dp.f [📝 src](fortran/int2p_dp.f)
**功能**: 在不同气压坐标层之间进行变量插值

**详细说明**:
- 包含一个子程序:
  - DINT2P: 气压坐标双精度插值程序
- 输入参数:
  - PPIN,XXIN: 输入气压和变量数组(NPIN)
  - PPOUT: 目标输出气压数组(NPOUT)
  - LINLOG: 插值类型(1=线性,0=对数线性,负值=允许外推)
  - XMSG: 缺失值标识
- 输出参数:
  - XXOUT: 插值后的变量数组
  - IER: 错误状态码
- 功能特点:
  - 支持递增和递减气压序列自动识别和重排序
  - 线性插值: $\text{SLOPE} = \frac{X_1-X_2}{P_1-P_2}$
  - 对数插值: $\text{SLOPE} = \frac{X_1-X_2}{\ln(P_1)-\ln(P_2)}$
  - 外推选项: 使用端点斜率进行边界外推
  - 自动处理缺失值和精确匹配的压力层
  - 适用于大气模式垂直坐标转换

## 72. julGreg.f [📝 src](fortran/julGreg.f)
**功能**: 儒略日数和格里高利历日期之间的相互转换

**详细说明**:
- 包含四个函数:
  - GREG2JULI: 格里高利历转整型儒略日数
  - GREG2JULD: 格里高利历转双精度儒略日数(含小时)
  - JULD2GREG: 双精度儒略日数转格里高利历
  - JULI2GREG: 整型儒略日数转格里高利历
- 输入参数:
  - GREG2JULI/JULD: YYYY,MM,DD[,HR](年,月,日[,时])
  - JULD2GREG/JULI2GREG: JULD/JULD(儒略日数)
- 输出参数:
  - 格里高利历: YYYY,MM,DD,HR
  - 儒略日数: INTEGER或REAL*8
- 功能特点:
  - 基于Fliegel & Van Flandern(1968)算法
  - 儒略日从公元前4713年1月1日起算
  - 儒略日以中午12:00 UT为起点
  - 支持小时精度的时间转换
  - 适用于历史和天文时间计算

## 73. kens_trend_dp.ncl.f [📝 src](fortran/kens_trend_dp.ncl.f)
**功能**: Kendall's $\tau$非参数趋势检验和Theil-Sen斜率估计

**详细说明**:
- 包含两个子程序:
  - KENSTSTNCL: NCL接口程序,处理工作数组分配
  - KENSTST: 核心Kendall趋势检验实现
- 输入参数:
  - XDATA: 时间序列数据(N)
  - N: 数据点数(必须$\geq$10)
  - EPS: 相等性判断阈值(默认1d-5)
- 输出参数:
  - S: Kendall's S统计量
  - Z: 标准化Z值
  - PROB: 显著性概率
  - SLOPE: Theil-Sen斜率数组(用于中位数斜率估计)
- 功能特点:
  - S统计量: $S = \sum_{j>i} \text{sign}(x_j - x_i)$
  - Z统计量: $Z = \frac{S \pm 1}{\sqrt{\text{Var}(S)}}$，方差调整并列值
  - 显著性: $\text{PROB} = \text{erf}\left(\frac{|Z|}{\sqrt{2}}\right)$
  - Theil-Sen斜率: 所有点对斜率的中位数
  - 适用于非正态分布和含异常值的趋势分析

## 74. kernel_density.f [📝 src](fortran/kernel_density.f)
**功能**: 高斯核密度估计和插件带宽选择算法

**详细说明**:
- 包含两个子程序:
  - KERDENI: NCL接口程序,处理缺失值和数据预处理
  - PLUGIN: 核心插件带宽选择和密度估计算法
- 输入参数:
  - X: 输入数据数组(N)
  - Z: 估计网格点数组(M)
  - XMSG: 缺失值标识
- 输出参数:
  - F: 各网格点的密度估计值(M)
  - H: 优化带宽值
- 功能特点:
  - 高斯核函数: K(u) = exp(-u$^2$/2)/$\sqrt{2\pi}$
  - 插件带宽选择: 迭代优化求解最优带宽
  - 密度估计: $f(x) = \frac{1}{nh}\sum K\left(\frac{x-X_i}{h}\right)$
  - IQR初始带宽: $h_0 = 0.92 \times \text{IQR} \times n^{-1/7}$
  - 5次迭代优化获得最终带宽
  - 适用于连续分布的非参数密度估计

## 75. kmeans_kmns_as136.f [📝 src](fortran/kmeans_kmns_as136.f)
**功能**: K-means聚类算法(AS-136标准实现)

**详细说明**:
- 包含四个子程序:
  - KMNS136: NCL接口程序,处理初始化和聚类中心设置
  - KMNS: K-means算法主程序
  - OPTRA: 最优传递阶段实现
  - QTRAN: 快速传递阶段实现
- 输入参数:
  - DAT: 数据矩阵(M×N),M个点,N维特征
  - K: 聚类数目
  - ITER: 最大迭代次数
  - ISEED: 初始化种子(1=顺序,2=间隔采样)
- 输出参数:
  - IC1: 各点的聚类分配(M)
  - CLCNTR: 最终聚类中心(N×K)
  - WSS: 各聚类内平方和(K)
  - IER: 错误码(0=成功,1=空聚类,2=未收敛,3=参数错误)
- 功能特点:
  - 基于Hartigan & Wong(1979)AS-136算法
  - 两阶段迭代: 最优传递(全局优化)+快速传递(局部优化)
  - 目标函数: 最小化类内平方和$\sum\sum||x_i - c_k||^2$
  - 自动检测空聚类和收敛性
  - 适用于数值数据的无监督聚类分析


## 76. kron_square.f [📝 src](fortran/kron_square.f)
**功能**: 计算两个方形矩阵的Kronecker积

**详细说明**:
- 包含一个子程序:
  - KRONSQM: 方形矩阵Kronecker积计算程序
- 输入参数:
  - N: 输入矩阵维度(N×N)
  - N2: 输出矩阵维度(N$^2$×N$^2$)
  - A,B: 两个输入方形矩阵(N×N)
- 输出参数:
  - C: Kronecker积矩阵(N$^2$×N$^2$)
- 功能特点:
  - Kronecker积定义: $C = A \otimes B$
  - 元素计算: $C(\text{row},\text{col}) = A(i,j) \times B(k,l)$
  - 索引映射: $\text{row} = n \times (i-1) + k$, $\text{col} = n \times (j-1) + l$
  - 输出矩阵维度为输入矩阵维度的平方
  - 适用于线性代数、矩阵分析、系统理论

## 77. lclvl.f [📝 src](fortran/lclvl.f)
**功能**: 计算气块抬升凝结高度(LCL)对应的压力值

**详细说明**:
- 包含一个子程序:
  - DLCLPRS: 抬升凝结高度压力计算程序
- 输入参数:
  - P: 初始压力(mb/hPa)
  - TK: 初始温度(K)
  - TDK: 露点温度(K)
- 输出参数:
  - PLCL: 抬升凝结高度对应的压力(mb/hPa)
- 功能特点:
  - 基于Stipanuk(1973)迭代算法
  - 饱和水汽压: $\text{ES} = \text{ES0} \times \exp\left(\frac{17.67 \times \text{TD}}{\text{TD}+243.5}\right)$
  - 饱和混合比: $\text{WS} = \frac{622 \times \text{ES}}{P - \text{ES}}$  
  - 位温: $\text{PTK} = \text{TK} \times \left(\frac{1000}{P}\right)^{0.286}$
  - 10次迭代求解干绝热线和等混合比线交点
  - 适用于大气物理、对流分析、天气预报

## 78. localMinMax.f [📝 src](fortran/localMinMax.f)
**功能**: 在二维网格中检测局地最小值和最大值

**详细说明**:
- 包含两个子程序:
  - DLOCALMN: 局地最小值检测程序
  - DLOCALMX: 局地最大值检测程序
- 输入参数:
  - X: 2D数据网格(NI×NJ)
  - XMSG: 缺失值标识
  - LWRAP: 边界环绕标志(1=经度环绕,0=不环绕)
  - DELTA: 极值判断阈值(通常=0.0)
- 输出参数:
  - XI,YJ: 极值点的i,j坐标(调整为NCL下标)
  - LOCMIN/LOCMAX: 极值点的数值
  - NMIN/NMAX: 检测到的极值点数量
- 功能特点:
  - 九点模板检测: 中心点与周围8个邻居比较
  - 最小值条件: 所有邻点 > 中心值+DELTA
  - 最大值条件: 所有邻点 < 中心值-DELTA
  - 支持经度环绕边界处理
  - 跳过包含缺失值的区域
  - 适用于气象场极值分析、模式识别

## 79. lspoly.f [📝 src](fortran/lspoly.f)
**功能**: 加权最小二乘多项式拟合算法

**详细说明**:
- 包含两个子程序:
  - DLSPOLY: NCL接口程序,处理参数检查和工作数组
  - DLSPLY2: 核心最小二乘拟合算法实现
- 输入参数:
  - X,Y: 数据点坐标(M个点)
  - WGT: 权重数组(M)
  - N: 多项式系数个数(度数=N-1)
  - M: 数据点数(必须M>N)
- 输出参数:
  - COEF: 多项式系数(N),COEF(1)为常数项
  - IER: 错误码(0=成功,1=M<N,2=N>5警告,3=矩阵奇异)
- 功能特点:
  - 正规方程法: $\mathbf{A}^T \mathbf{W} \mathbf{A} \times \text{COEF} = \mathbf{A}^T \mathbf{W} \mathbf{Y}$
  - 高斯消去法求解: 全选主元策略
  - 支持0-4度多项式(N$\leq$5)
  - 加权拟合: 权重为0表示排除该点
  - 数值稳定性: 高度多项式建议用其他方法
  - 适用于趋势拟合、曲线回归、数据平滑

## 80. mixhum_ptd.f [📝 src](fortran/mixhum_ptd.f)
**功能**: 根据压力、温度、露点温度计算水汽混合比或比湿

**详细说明**:
- 包含一个子程序:
  - DWMRQ: 水汽混合比/比湿计算程序
- 输入参数:
  - P: 压力数组(Pa),转换为mb/hPa
  - TD: 露点温度数组(K),转换为°C
  - PMSG,TDMSG: 压力和露点温度缺失值标识
  - ISWIT: 计算选项(1=混合比kg/kg, 2=比湿kg/kg, 负值=返回g/kg)
- 输出参数:
  - WMR: 水汽混合比或比湿数组
  - WMRMSG: 输出缺失值标识
- 功能特点:
  - 调用DWMRSKEWT函数计算基本混合比
  - 混合比(kg/kg): $\text{WMR} = \text{DWMRSKEWT}(P_{\text{mb}},TD_C) \times 0.001$
  - 比湿转换: $\text{SH} = \frac{\text{WMR}}{1 + \text{WMR}}$
  - 单位转换: 负ISWIT返回g/kg单位
  - 自动处理缺失值
  - 适用于大气湿度分析、水汽输送计算  
## 81. mixhum_ptrh.f [📝 src](fortran/mixhum_ptrh.f)
**功能**: 根据压力、温度、相对湿度计算水汽混合比或比湿

**详细说明**:
- 包含一个子程序:
  - DMIXHUM1: 水汽混合比/比湿计算程序
- 输入参数:
  - P: 压力(hPa或mb)
  - TK: 温度(K)
  - RH: 相对湿度(百分比)
  - ISWIT: 计算选项(1=混合比, 2=比湿, 负值=kg/kg单位)
- 输出参数:
  - QW: 水汽混合比或比湿
- 功能特点:
  - 饱和水汽压: $\text{EST} = 6.11 \times \exp\left(\frac{17.269 \times (\text{TK}-273.15)}{\text{TK}-35.86}\right)$
  - 混合比公式: $\text{QW} = \frac{0.622 \times \text{EST}}{P - 0.378 \times \text{EST}} \times \frac{\text{RH}}{100}$
  - 比湿转换: $\text{SH} = \frac{\text{QW}}{1 + \text{QW}}$
  - 单位控制: 负ISWIT返回kg/kg, 正值返回g/kg
  - 适用于大气湿度计算、露点分析

## 82. mlegev_memory.f [📝 src](fortran/mlegev_memory.f)
**功能**: 使用AS-215算法进行广义极值分布参数的最大似然估计

**详细说明**:
- 包含两个子程序:
  - DMLEGEVI: NCL接口程序,处理缺失值和参数初值设定
  - DMLEGEV: AS-215核心算法实现
- 输入参数:
  - X: 数据数组(去除缺失值后$\geq$10个点)
  - XMSG: 缺失值标识
  - MONIT: 监控标志(0=无监控)
- 输出参数:
  - VALS(6): [位置参数$\xi$, 尺度参数$\sigma$, 形状参数$\kappa$, 标准误差×3]
  - IERR: 错误码(0=成功, 1=N<2, 2=未收敛, 3=超限, 4=步长限制)
- 功能特点:
  - 基于Hosking(1985)AS-215算法
  - Newton-Raphson和最陡上升混合迭代
  - 自动参数初值估计: $\sigma_0=\text{std}\times\sqrt{6/\pi}$, $\xi_0=\text{mean}-\gamma\times\sigma_0$
  - 最大30次迭代,50次似然评估
  - 适用于极值统计、风险分析、水文频率分析

## 83. obsp1_mult_time_dp.f [📝 src](fortran/obsp1_mult_time_dp.f)
**功能**: 基于迭代改进的客观分析插值算法

**详细说明**:
- 包含两个子程序:
  - DOBJANLX: 驱动程序,处理数据预处理和排序
  - DOBJANL: 核心客观分析算法实现
- 输入参数:
  - PLON,PLAT,PVAL: 观测点经纬度和数值(NTIM×NPTS)
  - GLAT,GLON: 输出网格经纬度(NLAT×MLON)
  - RSCAN: 扫描半径序列(度),递减排列
  - SMWGT: 混合权重(NSCAN)
  - XMSG,PMSG: 缺失值标识
- 输出参数:
  - GRID: 插值网格(MLON×NLAT×NTIM)
  - IER: 错误码
- 功能特点:
  - 多尺度分析: 大→小半径逐级细化
  - 高斯权重: wgt = exp(-4×(距离/半径)$^2$)
  - 距离计算: 大圆距离,精确球面几何
  - 迭代混合: grid = smwgt×新值 + (1-smwgt)×旧值
  - 适用于气象客观分析、台站数据格点化

## 84. ocean.f [📝 src](fortran/ocean.f)
**功能**: 海洋数据的区域加权平滑和混合层深度计算

**详细说明**:
- 包含两个子程序:
  - WGT_AREA_SMOOTH: 区域加权平滑算法
  - MIXED_LAYER_DEPTH: 混合层深度计算
- WGT_AREA_SMOOTH输入参数:
  - FIELD: 输入场(NX×NY×NO)
  - AREA: 面积权重(NX×NY)
  - ICYCLIC: 经度循环边界标志
  - FILL_VALUE: 缺失值
- MIXED_LAYER_DEPTH输入参数:
  - FIELD: 密度场(NX×NY×NZ)
  - KMT: 海底地形(NX×NY)
  - HT: 海表高度(NX×NY)
  - DEPTH: 深度数组(NZ)
  - OFFSET: 密度判别阈值
- 输出参数:
  - FIELD_RET: 平滑后场或混合层深度
- 功能特点:
  - 五点加权平滑: 中心+东西南北邻点
  - 面积权重归一化: 权重和=1
  - 混合层定义: 表层密度+offset的深度位置
  - 线性插值求深度: depth = d₁ + ($\rho$-$\rho$₁)×(d₂-d₁)/($\rho$₂-$\rho$₁)
  - 适用于海洋环流分析、混合层研究

## 85. omcalc_ccm.f [📝 src](fortran/omcalc_ccm.f)
**功能**: CCM气候模式中垂直压力速度(omega)的诊断计算

**详细说明**:
- 包含一个子程序:
  - OMCALCCCM: 基于CAM3.0修改的omega计算程序
- 输入参数:
  - U,V: 纬向和经向风速(ILON×JLAT×KLEV)
  - D: 散度场(ILON×JLAT×KLEV)
  - DPSL,DPSM: 地面压力梯度经纬度分量
  - PMID,PDEL: 层中压力和层厚(ILON×JLAT×KLEV)
  - PSFC: 地面压力(ILON×JLAT)
  - HYBD,HYBM: 混合坐标系数(KLEV)
  - NPRLEV: 第一个纯压力层
- 输出参数:
  - OMEGA: 垂直压力速度(ILON×JLAT×KLEV)
- 功能特点:
  - 先计算$\omega$/p,后缩放为$\omega$
  - 静力学矩阵方程: HKK×$\omega$/p + HLK×$\sum$散度×Δp
  - 混合坐标贡献: $\vec{v} \cdot \nabla p_s$项的处理
  - 层间积分: 从顶层到底层逐层计算
  - 最终缩放: $\omega$ = ($\omega$/p)×p_mid
  - 适用于大气动力学诊断、垂直运动分析  
## 86. p2hyo.f [📝 src](fortran/p2hyo.f)
**功能**: 将常压面数据插值到混合坐标层

**详细说明**:
- 包含两个子程序:
  - P2HYO: 主插值接口程序，处理参数检查和数组检验
  - P2HYB: 核心插值算法实现
- 输入参数:
  - PI: 输入压力层数组(KLEVI)，必须单调
  - XI: 输入变量(MLON×NLAT×KLEVI)
  - PSFC: 地面压力(MLON×NLAT)，单位Pa
  - P0: 基准压力，HYAO系数的参考压力
  - HYAO,HYBO: 输出混合坐标系数(KLEVO)
  - KFLAG: 外推选项(0-4，控制边界外处理方式)
- 输出参数:
  - XO: 插值后变量(MLON×NLAT×KLEVO)
  - IFLAG: 缺失值存在标志
  - IER: 错误码(0=成功，其他=压力数组非单调)
- 功能特点:
  - 混合坐标压力: p(k) = hyao(k)×p0 + hybo(k)×psfc
  - 对数线性插值: 在ln(p)空间内线性插值
  - 5种外推模式: 无外推、最近值、单边外推、双边外推
  - 边界处理: DXDP梯度外推公式
  - 适用于气候模式垂直坐标转换

## 87. paleo_coasts.f [📝 src](fortran/paleo_coasts.f)
**功能**: 从掩膜数据生成EZMAP古地理海岸线数据库

**详细说明**:
- 包含多个子程序:
  - PALEOOUTLINE: 主程序，协调掩膜转换和文件生成
  - SVBLED: 核心边界线追踪算法
  - WTEMLR: EZMAP线数据文件写入
  - WRNUMB,WRCHAR: 数字和字符格式化输出
- 输入参数:
  - MASK: 陆海掩膜数组(NLON×NLAT)，MSKVAL=陆地
  - LAT,LON: 纬度和经度数组(NLAT,NLON)
  - ZDAT: 临时数组(IM×JM)，IM=2×NLON+1，JM=2×NLAT+1
  - FLINES,FNAMES: 输出的线数据和名称数据文件名
  - MSKVAL: 陆地掩膜值
- 输出参数:
  - 生成的EZMAP数据库文件
- 功能特点:
  - 边界追踪: 8方向矢量搜索算法
  - 坐标变换: 掩膜网格到地理坐标的线性映射
  - 文件格式: 标准EZMAP二进制格式
  - 区域命名: 左侧陆地("Land")，右侧海洋("Water")
  - 适用于古气候研究、地质历史重建

## 88-92. PASSB系列 - FFT后向传递函数集
**功能**: FFTPACK中FFT后向变换(逆变换)的传递函数

**详细说明**:
- PASSB: 通用后向传递函数，处理任意质因子
- 输入参数:
  - NAC,IDO,IP,L1,IDL1: FFT算法控制参数
  - CC,C1,C2,CH,CH2: 输入输出数组
  - WA: 预计算的三角函数表
- 输出参数:
  - 变换后的复数数组
- 功能特点:
  - 混合基算法: 处理2,3,4,5等质因子分解
  - 双缓冲技术: CC↔CH数组交替使用
  - 三角函数优化: 预存储避免重复计算
  - 对称性利用: IPPH=(IP+1)/2半长度计算
  - 适用于复数到时域的FFT逆变换

## 93-97. PASSF系列 - FFT前向传递函数集  
**功能**: FFTPACK中FFT前向变换的传递函数

**详细说明**:
- PASSF: 通用前向传递函数，处理任意质因子
- 输入参数:
  - NAC,IDO,IP,L1,IDL1: FFT算法控制参数
  - CC,C1,C2,CH,CH2: 输入输出数组
  - WA: 预计算的三角函数表
- 输出参数:
  - 变换后的复数数组
- 功能特点:
  - 与PASSB算法对称: 符号差异在关键运算处
  - 前向变换: 时域到频域转换
  - 高效分解: 大长度FFT分解为小长度FFT乘积
  - 位反转: 自然顺序到位反转顺序转换
  - 适用于时域到复数频域的FFT变换

## 98. patternCor.f [📝 src](fortran/patternCor.f)
**功能**: 计算两个二维场之间的模式(异常)相关系数

**详细说明**:
- 包含两个子程序:
  - PATCOR1: 使用一维权重的模式相关计算
  - PATCOR2: 使用二维权重的模式相关计算
- 输入参数:
  - X,Y: 两个数据场(MLON×NLAT)
  - W: 权重数组(NLAT或MLON×NLAT)
  - XMSG,YMSG: 缺失值标识
- 输出参数:
  - R: 模式相关系数
  - IER: 错误码(0=成功,1=全缺失,2=常数场)
- 功能特点:
  - 加权平均: $\bar{x} = \frac{\sum(w \times x)}{\sum w}$, $\bar{y} = \frac{\sum(w \times y)}{\sum w}$
  - 异常场: $x' = x-\bar{x}$, $y' = y-\bar{y}$
  - 相关系数: $r = \frac{\sum(w \times x' \times y')}{\sqrt{\sum(w \times x'^2) \times \sum(w \times y'^2)}}$
  - 权重支持: 纬度权重(cos φ)或面积权重
  - 适用于气候模式评估、场相似性分析

## 99. patternCor3.f [📝 src](fortran/patternCor3.f)
**功能**: 计算三维数据场之间的模式相关系数

**详细说明**:
- 扩展自patternCor.f的三维版本
- 功能特点:
  - 支持三维权重分布
  - 垂直层间相关分析
  - 时空模式识别
  - 多层大气/海洋数据相关性
  - 适用于3D气候场对比、模式验证

## 100. phybrid_ccm.f [📝 src](fortran/phybrid_ccm.f)
**功能**: CCM气候模式混合坐标系统的压力场计算

**详细说明**:
- 基于CCM/CAM混合垂直坐标系统
- 功能特点:
  - 混合坐标公式: $p(i,j,k) = \text{hyam}(k) \times p_0 + \text{hybm}(k) \times p_s(i,j)$
  - 层边界压力: $p_{\text{interface}} = \text{hyai}(k) \times p_0 + \text{hybi}(k) \times p_s(i,j)$
  - 层厚计算: $\Delta p = p_{\text{interface}}(k+1) - p_{\text{interface}}(k)$
  - 地形跟随: 近地面$\sigma$坐标，高空压力坐标
  - 适用于大气模式诊断、垂直坐标处理  
## 101. plotfmt_rddata.f [📝 src](fortran/plotfmt_rddata.f)
**功能**: 绘图格式数据读取接口

**详细说明**:
- 包含一个子程序:
  - PLOTFMT_RDDATA: 从固定文件单元读取二维数据数组
- 输入参数:
  - NX,NY: 数据数组维度
  - FUNIT: 固定为10的文件单元号
- 输出参数:
  - SLAB: 读取的二维数据数组(NX×NY)
  - ISTATUS: 状态码(0=成功,1=读取错误)
- 功能特点:
  - 使用Fortran无格式二进制读取
  - 简单的错误处理机制
  - 固定文件单元号设计
  - 适用于绘图程序的数据输入接口

## 102. potmp_dpth2pres.f [📝 src](fortran/potmp_dpth2pres.f)
**功能**: 海洋学位温计算和深度到压力转换

**详细说明**:
- 包含两个子程序:
  - DPOTMP: 海水位温计算
  - DPTH2PRES: 深度到压力坐标转换
- DPOTMP输入参数:
  - PRESS: 压力(分巴)
  - TEMP: 温度(摄氏度)
  - S: 盐度(PSS-78标准)
  - RP: 参考压力(分巴)
- DPTH2PRES输入参数:
  - DEPTH: 深度(米)
  - HAS_MSG_DEPTH: 缺失值标志
  - DEPTH_MSG: 深度缺失值代码
- 输出参数:
  - POTEMP: 位温(摄氏度)
  - PRESSURE: 压力(巴)
- 功能特点:
  - 基于Fofonoff(1976)位温计算公式
  - 使用Levitus(1994)全球平均温盐数据
  - 4阶Runge-Kutta积分方法
  - 静力平衡原理: $dp = -\rho g dz$
  - 压力公式: $p = 0.059808 \times (e^{-0.025z}-1) + 0.100766z + 2.28405 \times 10^{-7} \times z^2$
  - 适用于海洋学物理海洋分析

## 103. prcwat.f [📝 src](fortran/prcwat.f)
**功能**: 可降水量(柱总水汽)计算

**详细说明**:
- 包含一个子程序:
  - DPRCWATDP: 从比湿和层厚计算可降水量
- 输入参数:
  - Q: 比湿数组(kg/kg，KLVL层)
  - DP: 层厚数组(Pa，KLVL层)
  - KLVL: 垂直层数
  - QMSG,DPMSG: 比湿和层厚的缺失值代码
- 输出参数:
  - PRCWAT: 可降水量(kg/m$^2$或mm)
- 功能特点:
  - 积分公式: PRCWAT = (1/g)$\sum$[q(k)×|dp(k)|]
  - 使用标准重力加速度: g = 9.81 m/s$^2$
  - 自动处理缺失值层
  - 单位换算: kg/m$^2$ = mm水柱高度
  - 适用于大气水汽分析、降水潜力评估

## 104. pres_hybrid_jra55_dp.f [📝 src](fortran/pres_hybrid_jra55_dp.f)
**功能**: JRA55再分析混合坐标压力计算

**详细说明**:
- 包含一个子程序:
  - DPHYBRIDJRA55: JRA55专用混合坐标压力计算
- 输入参数:
  - HYA,HYB: JRA55界面混合坐标系数(KLEVI层)
  - PSFC: 地面压力(MLON×NLAT，Pa)
  - MLON,NLAT,KLEV,KLEVI: 网格和层数维度
  - PI: 界面压力工作数组
  - PMSG: 缺失值代码
- 输出参数:
  - PHY: 层中压力(MLON×NLAT×KLEV，Pa)
  - PI: 界面压力(MLON×NLAT×KLEVI，Pa)
- 功能特点:
  - 界面压力: $p_i(k) = \text{hya}(k) + \text{hyb}(k) \times p_{\text{sfc}}$
  - 层中压力对数加权公式: $p = \exp\left[\frac{1}{\Delta p} \times (p_1 \ln p_1 - p_2 \ln p_2) - 1\right]$
  - $\Delta p = p_{\text{interface}}(k) - p_{\text{interface}}(k+1)$
  - 底层强制设为10Pa
  - 基于JRA55再分析系统文档
  - 适用于JRA55数据的垂直坐标处理

## 105. preshybrid.f [📝 src](fortran/preshybrid.f)
**功能**: 混合坐标压力计算(通用版本)

**详细说明**:
- 包含两个子程序:
  - PRESHYBRID: 混合坐标层中压力计算
  - DPRESHYBRID: 混合坐标层厚计算
- PRESHYBRID输入参数:
  - $P_0$: 基准压力 (Pa)
  - $P_s$: 地面压力 (Pa)
  - HYA, HYB: 混合坐标系数 (KLVL 层)
- DPRESHYBRID输入参数:
  - $P_0$, $P_s$: 基准压力和地面压力 (Pa)
  - HYAI, HYBI: 界面混合坐标系数 (KLVL 层)
- 输出参数:
  - PHY: 层中压力 (KLVL, Pa)
  - DPHY: 层厚 (KLVL-1, Pa)
- 功能特点:
  - 标准混合坐标公式: $p(k) = \text{hya}(k) \times P_0 + \text{hyb}(k) \times P_s$
  - 层厚计算: $\Delta p(k) = \left| P_0 \times [\text{hyai}(k+1) - \text{hyai}(k)] + [\text{hyai}(k+1) - \text{hybi}(k)] \times P_s \right|$
  - 支持 CCM2/CCM3 混合坐标系统
  - 从模式顶到地面排序
  - 适用于大气模式垂直坐标转换

## 106. prneof_dp.f [📝 src](fortran/prneof_dp.f)
**功能**: 主成分EOF分析(双精度LAPACK实现)

**详细说明**:
- 包含两个子程序:
  - DDRVEOF: 主EOF分析驱动程序
  - DNCLDRV: NCL接口驱动程序，带缺失值阈值控制
- 输入参数:
  - X: 时空数据矩阵(NOBS×MSTA)
  - NOBS,MSTA: 观测数和站点数
  - XMSG: 缺失值代码
  - NEVAL: 需要计算的特征值/特征向量个数
  - JOPT: 矩阵类型(0=协方差矩阵,1=相关矩阵)
  - PCRIT: 有效数据百分比阈值(DNCLDRV)
- 输出参数:
  - EVAL: 特征值(降序排列)
  - EVEC: 特征向量(空间模态)
  - PCVAR: 各模态方差贡献百分比
  - TRACE: 协方差/相关矩阵的迹
- 功能特点:
  - 基于LAPACK DSPEVX高精度特征值求解
  - 支持协方差和相关矩阵两种分析
  - 自动处理缺失值和数据质量控制
  - 特征值按方差贡献降序排列
  - 可计算主成分时间系数(需额外计算)
  - 适用于气候模态分析、降维、模式识别

## 107. prneofTranspose.f [📝 src](fortran/prneofTranspose.f)
**功能**: EOF分析矩阵转置优化版本

**详细说明**:
- 包含两个子程序:
  - TDRVPRC: 转置EOF分析主接口
  - XRVEOFT: 转置矩阵EOF核心算法
- 输入参数:
  - X: 原始数据矩阵(NROW×NCOL)
  - NROBS,NCSTA: 观测数和站点数
  - XMSG: 缺失值代码
  - NEVAL: 特征值个数
  - PCRIT: 数据完整性阈值百分比
  - JOPT: 标准化选项(0=协方差,1=相关)
- 输出参数:
  - EVAL: 特征值
  - EVEC: 重构的空间特征向量
  - PCVAR: 方差贡献百分比
  - TRACE: 协方差矩阵迹
- 功能特点:
  - 基于数据矩阵转置提高计算效率
  - 当NOBS<<NCSTA时显著加速计算
  - 利用时空对偶性: EOF_space = f(PC_time, Data_transpose)
  - 自动数据标准化和异常值处理
  - 质量控制: 剔除有效数据<PCRIT的站点
  - 自动分配工作数组，无需预设内存
  - 适用于大规模网格数据EOF分析

## 108. psigma.f [📝 src](fortran/psigma.f)
**功能**: $\sigma$坐标系统压力计算

**详细说明**:
- 包含一个子程序:
  - DPSIGMA: $\sigma$坐标到压力坐标转换
- 输入参数:
  - SIGMA: $\sigma$坐标系数(KLEV层)
  - PSFC: 地面压力(MLON×NLAT，Pa)
  - MLON,NLAT,KLEV: 经度、纬度、层数维度
- 输出参数:
  - PSIG: $\sigma$层压力(MLON×NLAT×KLEV，Pa)
- 功能特点:
  - 简单$\sigma$坐标公式: $p(i,j,k) = \sigma(k) \times p_{\text{sfc}}(i,j)$
  - $\sigma$坐标定义: $\sigma = p/p_{\text{sfc}}$，范围[0,1]
  - $\sigma$=0对应模式顶，$\sigma$=1对应地面
  - 地形跟随坐标系统
  - 保持垂直分辨率随地形变化
  - 适用于地形复杂区域的大气模式
  - 比混合坐标系统更简单直接

## 109. random_dp.f [📝 src](fortran/random_dp.f)
**功能**: 双精度随机数生成库(完整统计分布)

**详细说明**:
- 包含多个随机数生成函数:
  - DGENCHI: 卡方分布随机数生成
  - DGENGAM: 伽马分布随机数生成  
  - DGENNOR: 正态分布随机数生成
  - DGENUNF: 均匀分布随机数生成
  - IGNLGI: 大整数均匀随机数生成
  - IGNUIN: 指定范围整数随机数
- 基础随机数引擎:
  - DRANFD: 基础[0,1)均匀随机数生成器
  - 基于L'Ecuyer & Cote(1991)算法
  - 支持32个独立随机数流
  - 周期长度约2.3×10^18
- 管理函数:
  - SETALL/SETSD: 设置随机数种子
  - GETSD: 获取当前种子状态
  - SETCGN/DGETCGN: 切换随机数生成器
- 功能特点:
  - 高质量L'Ecuyer线性同余生成器
  - 完整的概率分布函数库
  - 支持多流并行随机数生成
  - 通过统计验证的算法实现
  - 适用于Monte Carlo模拟、统计抽样

## 110. rdsstoi.f [📝 src](fortran/rdsstoi.f)
**功能**: NOAA最优插值海表温度数据读取

**详细说明**:
- 包含两个子程序:
  - RDSSTOI: 读取SSTOI格式的SST数据
  - WRSSTOI: 写入高斯网格SST数据
- RDSSTOI输入参数:
  - NYRSTRT,NYRLAST: 数据年份范围
  - MLON,NLAT: 网格维度(固定360×180)
- RDSSTOI输出参数:
  - SST: 海表温度数组(360×180，°C)
  - INFO: 日期和状态信息(9个整数)
- 数据格式特征:
  - 1°×1°全球网格(360经度×180纬度)
  - 地理范围: 179.5°W-179.5°E, 89.5°S-89.5°N
  - 支持周、月和气候态数据
  - 原始数据为°C×100的整数格式
  - 自动转换为浮点温度值
- 功能特点:
  - 顺序读取formatted文件'SSTOI'
  - 每个网格包含8项头信息
  - 头信息格式: (起始日期,结束日期,天数,索引)
  - 数据质量标识: $\leq$-1.78°C表示海冰
  - 需配合陆海掩膜使用
  - 适用于海洋学分析、气候研究

## 111. regcoef_dp.f [📝 src](fortran/regcoef_dp.f)
**功能**: 线性回归系数计算和统计检验

**详细说明**:
- 包含一个子程序:
  - DREGCOEF: 完整的线性回归分析
- 输入参数:
  - X,Y: 输入数据向量(NPTS个点)
  - NPTS: 数据点总数
  - XMSG,YMSG: X和Y的缺失值代码
- 输出参数:
  - RCOEF: 回归系数(斜率)
  - TVAL: t统计量(检验回归系数显著性)
  - NPTXY: 实际使用的数据点数
  - XAVE,YAVE: X和Y的均值
  - RSTD: 回归系数标准误差
  - YINT: y截距
  - IER: 错误代码
- 功能特点:
  - 基于Brownlee(1965)统计理论
  - 回归方程: $Y = \text{RCOEF} \times X + \text{YINT}$
  - t统计量: $t = \frac{\text{RCOEF}-0}{\text{RSTD}}$，用于显著性检验
  - 自由度: $df = n-2$
  - 自动处理缺失值和边界条件
  - 包含完整的回归诊断统计量
  - 适用于趋势分析、相关性研究

## 112. relhum_dp.f [📝 src](fortran/relhum_dp.f)
**功能**: 相对湿度计算(双精度，基于CCM算法)

**详细说明**:
- 包含一个函数:
  - DRELHUM: 从温度、混合比、压力计算相对湿度
- 输入参数:
  - T: 温度(K)
  - W: 混合比(kg/kg)
  - P: 压力(Pa)
- 输出参数:
  - DRELHUM: 相对湿度(%)
- 功能特点:
  - 内置饱和水汽压查找表(-100°C到+102°C)
  - 温度范围限制: 173.16K到375.16K
  - 相对湿度公式: $\text{RH} = \frac{W \times (P - 0.378 \times \text{Es})}{0.622 \times \text{Es}} \times 100\%$
  - 线性插值计算饱和水汽压
  - 下限约束: RH$\geq$0.0001%
  - 基于CCM气候模式物理过程
  - 适用于大气湿度分析、模式诊断

## 113. relhum_ice.f [📝 src](fortran/relhum_ice.f)
**功能**: 冰面相对湿度计算(改进的Magnus公式)

**详细说明**:
- 包含一个子程序:
  - DRELHUMI: 相对于冰面的相对湿度计算
- 输入参数:
  - P: 压力(Pa，内部转换为hPa)
  - TK: 温度(K)
  - QW: 混合比(kg/kg)
- 输出参数:
  - RH: 相对湿度(%)
- 功能特点:
  - 基于Alduchov & Eskridge改进Magnus公式
  - 适用温度范围: +50°C到-80°C
  - 最大相对误差: 0.337%-0.823%
  - 冰面饱和水汽压: $\text{Es} = 6.1128 \times \exp\left(\frac{22.571 \times (T-273.15)}{T-0.44}\right)$
  - 饱和混合比: $\text{Qs} = \frac{0.622 \times \text{Es}}{P - 0.378 \times \text{Es}}$
  - 相对湿度: $\text{RH} = \frac{100 \times \text{QW}}{\text{Qs}}$
  - 适用于低温条件、高空大气、极地气象

## 114. relhum_water.f [📝 src](fortran/relhum_water.f)
**功能**: 水面相对湿度计算(与mixhum_ptrh可逆计算)

**详细说明**:
- 包含一个子程序:
  - DRELHUMW: 相对于液态水的相对湿度计算
- 输入参数:
  - P: 压力(Pa，内部转换为hPa)
  - TK: 温度(K)  
  - QW: 混合比(kg/kg)
- 输出参数:
  - RH: 相对湿度(%)
- 功能特点:
  - 使用与mixhum_ptrh相同的常数，确保可逆计算
  - 水面饱和水汽压: $\text{Es} = 6.11 \times \exp\left(\frac{17.269 \times (T-273.15)}{T-35.86}\right)$
  - 饱和混合比: $\text{Qs} = \frac{0.622 \times \text{Es}}{P - 0.378 \times \text{Es}}$
  - 相对湿度: $\text{RH} = \frac{100 \times \text{QW}}{\text{Qs}}$
  - 分子量比: $M_w/M_d = 0.622$
  - 适用于正常大气条件、海洋大气、湿润环境

## 115. remap.f [📝 src](fortran/remap.f)
**功能**: 预计算权重的网格重映射

**详细说明**:
- 包含一个子程序:
  - DPOPREMAP: 使用预计算权重进行网格数据重映射
- 输入参数:
  - SRC_ARRAY: 源网格数据(NSRC)
  - MAP_WTS: 重映射权重矩阵(NW×NLINK)
  - DST_ADD,SRC_ADD: 目标和源网格地址索引(NLINK)
  - NLINK: 插值链接数
  - NW: 权重数(通常为1)
  - NDST,NSRC: 目标和源网格点数
  - XMSG: 缺失值代码
- 输出参数:
  - DST_ARRAY: 重映射后的目标网格数据(NDST)
- 功能特点:
  - 基于预计算的稀疏权重矩阵
  - 加权插值: DST(i) = $\Sigma$[SRC(j)×WGT(i,j)]
  - 支持任意网格间的映射关系
  - 自动处理缺失值传播
  - 高效的稀疏矩阵存储格式
  - 适用于保守插值、双线性插值、距离加权等
  - 常用于气候模式数据后处理

## 116. rhombtri_dp.f [📝 src](fortran/rhombtri_dp.f)
**功能**: 球谐系数菱形和三角形截断(双精度)

**详细说明**:
- 包含三个子程序:
  - DRHOMBTRUNC: 菱形截断球谐系数
  - DTRITRUNC: 三角形截断接口程序
  - DTRUNC: 三角形截断核心算法
- DRHOMBTRUNC输入参数:
  - AR,AI: 球谐系数实部和虚部(M×N)
  - M,N: 纬向和总波数维度
  - R: 菱形截断参数(保持的总波数)
- DTRUNC输入参数:
  - A,B: 球谐系数数组(ID×NM)
  - MS: 截断参数
  - ID,NM: 数组维度
- 输出参数:
  - 截断后的球谐系数(超出截断范围的系数置零)
- 功能特点:
  - 菱形截断: 在波数空间形成菱形保留区域
  - 三角形截断: 防止乘积项产生混叠
  - 基于SPHEREPACK 2.0算法
  - 支持不同M$\geq$N和M<N的情况
  - 用于球面谱方法中的去混叠
  - 适用于全球数值模式的谱变换

## 117. rmStndx_dp.f [📝 src](fortran/rmStndx_dp.f)
**功能**: 时间序列均值、中位数移除和标准化

**详细说明**:
- 包含三个子程序:
  - DRMVMEAN: 移除序列均值
  - DRMVMED: 移除序列中位数
  - DXSTND: 标准化(零均值单位方差)
- DRMVMEAN输入参数:
  - X: 输入时间序列(NPTS，原地修改)
  - NPTS: 数据点数
  - XMSG: 缺失值代码
- DXSTND输入参数:
  - X: 输入时间序列(NPTS，原地修改)
  - IOPT: 标准差类型(0=样本标准差，1=总体标准差)
- 输出参数:
  - X: 处理后的序列
  - IER: 错误代码
- 功能特点:
  - 原地操作，节省内存
  - 自动跳过缺失值
  - 标准化: X' = (X-μ)/$\sigma$
  - 支持样本和总体两种标准差计算
  - 适用于时间序列预处理、异常值标准化

## 118. rmsd.f [📝 src](fortran/rmsd.f)
**功能**: 均方根差(RMSD)计算

**详细说明**:
- 包含一个子程序:
  - DRMSD: 计算两个数组间的均方根差
- 输入参数:
  - X,Y: 两个输入向量(NPTS)
  - NPTS: 数据点总数
  - XMSG,YMSG: X和Y的缺失值代码
- 输出参数:
  - XYRMSD: 均方根差
  - NPTUSE: 实际使用的数据点数
  - IER: 错误代码
- 功能特点:
  - RMSD公式: $\\sqrt{}$[$\Sigma$(X-Y)$^2$/N]
  - 自动处理缺失值
  - 只有两个数组对应位置都有有效值才参与计算
  - 返回实际使用的点数用于验证
  - 适用于模式验证、预报评估、数据对比

## 119. round.f [📝 src](fortran/round.f)
**功能**: 数值四舍五入到最近整数

**详细说明**:
- 包含一个子程序:
  - RNDNCL: 将数组元素四舍五入到最近整数
- 输入参数:
  - XIN: 输入浮点数数组(NPTS)
  - NPTS: 数组长度
  - ISMSG: 是否有缺失值(0=无，1=有)
  - XMSG: 缺失值代码
  - IOPT: 输出选项(3=整数输出时调整缺失值)
- 输出参数:
  - XOUT: 四舍五入后的数组(NPTS)
- 功能特点:
  - 使用ANINT内置函数进行四舍五入
  - 自动处理缺失值传播
  - IOPT=3时自动调整非整数缺失值
  - 保持原数组不变(非原地操作)
  - 适用于数据离散化、整数转换

## 120. runave_dp.f [📝 src](fortran/runave_dp.f)
**功能**: 滑动平均滤波(双精度)

**详细说明**:
- 包含两个子程序:
  - DRUNAVE: 滑动平均主接口
  - DRUNAVX77: 滑动平均核心算法
- 输入参数:
  - X: 输入时间序列(NPTS，原地修改)
  - NPTS: 数据点数
  - NAVE: 滑动窗口大小
  - KOPT: 边界处理选项
  - XMSG: 缺失值代码
  - WORK: 工作数组(LWORK)
- 边界处理选项:
  - KOPT<0: 循环边界条件
  - KOPT=0: 边界点设为缺失值
  - KOPT>0: 反射(对称)边界条件
- 输出参数:
  - X: 平滑后的时间序列
  - IER: 错误代码
- 功能特点:
  - 支持奇数和偶数窗口大小
  - 三种边界处理方式
  - 自动处理序列中的缺失值
  - 窗口内任何缺失值导致该点输出为缺失值
  - 适用于时间序列平滑、噪声抑制

## 121. runknt.f [📝 src](fortran/runknt.f)
**功能**: 连续序列计数统计

**详细说明**:
- 包含两个子程序:
  - NUMRUN2: 二维数组连续序列计数
  - NUMRUN1: 一维数组连续序列计数  
- 输入参数:
  - XIN: 输入整数数组(0/1值)
  - NX,NTIM: 数组维度
  - OPTKNT: 计数选项(0=重叠计数，1=唯一序列计数)
- 输出参数:
  - XKNT: 各长度连续序列的计数(NX×NTIM或NX)
- 功能特点:
  - OPTKNT=0: 计数所有长度为n的连续1序列(重叠)
  - OPTKNT=1: 计数独立的连续1序列(两端必须为0)
  - 自动填充边界为0以确保序列完整计数
  - 适用于干旱/湿润期统计、连续事件分析
  - 气候极端事件持续时间统计

## 122. rvdv.f [📝 src](fortran/rvdv.f)
**功能**: 相对涡度和散度有限差分计算

**详细说明**:
- 包含四个子程序:
  - DVRFIDF: 相对涡度有限差分计算
  - DDVFIDF: 散度有限差分计算
  - DLNEXTRP: 角点线性外推
  - DUSEAVE: 极点平均值处理
- 输入参数:
  - U,V: 风场分量(MLON×NLAT，m/s)
  - GLAT,GLON: 纬度和经度数组
  - IOPT: 边界条件选项
  - XMSG: 缺失值代码
- 输出参数:
  - RV: 相对涡度(s$^{-1}$)
  - DV: 散度(s$^{-1}$)
- 功能特点:
  - 相对涡度: $\zeta = \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y} + \frac{u}{a}\tan(\phi)$
  - 散度: $\nabla \cdot \vec{V} = \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} - \frac{v}{a}\tan(\phi)$  
  - 支持等间距经度和不等间距纬度
  - 三种边界处理选项：循环、非循环、混合
  - 极点使用平均值，角点线性外推
  - 适用于大气动力学分析

## 123. s2h.f [📝 src](fortran/s2h.f)
**功能**: $\sigma$坐标到混合坐标插值

**详细说明**:
- 包含两个子程序:
  - DH2SDRV: $\sigma$到混合坐标转换驱动程序
  - DS2HBD: $\sigma$坐标插值核心算法
- 输入参数:
  - DATI: $\sigma$坐标数据(NLVI)
  - HYA,HYB: 混合坐标系数(NLVO)
  - P0: 基准压力(Pa)
  - PSFC: 地面压力(Pa)
  - SIGI: $\sigma$坐标层(NLVI)
  - INTYP: 插值类型(1=线性，2=对数，3=双对数)
- 输出参数:
  - DATO: 混合坐标数据(NLVO)
  - SIGO: 转换的$\sigma$值(NLVO)
- 功能特点:
  - 混合坐标转$\sigma$: $\sigma(k) = \frac{\text{hya}(k) \times P_0}{P_{\text{sfc}}} + \text{hyb}(k)$
  - 三种插值方法：线性、对数、双对数
  - 自动处理模式顶、底部和内部层
  - 适用于大气模式垂直坐标变换

## 124. setrig.f [📝 src](fortran/setrig.f)
**功能**: Lin-Rood模式纬度和权重计算

**详细说明**:
- 包含三个子程序:
  - LINROOD: Lin-Rood纬度和权重主接口
  - LINROODWT: Lin-Rood权重计算
  - SETRIG: 三角函数网格设置
- 输入参数:
  - NLAT: 纬度点数
- 输出参数:
  - LAT: 纬度数组(度)
  - WEIGHT: 高斯权重数组
  - COSP,COSE,SINP,SINE: 三角函数数组
- 功能特点:
  - 均匀纬度网格: $\Delta\text{lat} = \frac{180°}{\text{nlat}-1}$
  - 权重计算: $w(j) = \sin(\phi_{j+1}) - \sin(\phi_j)$
  - 极点特殊处理: $w(1) = 1 + \sin(\phi_2)$, $w(\text{nlat}) = 1 - \sin(\phi_{\text{nlat}})$
  - 基于NASA Lin-Rood动力框架
  - 适用于有限体积大气模式

## 125. sg_tools.f [📝 src](fortran/sg_tools.f)
**功能**: 球面几何计算工具集

**详细说明**:
- 包含多个球面几何函数:
  - GCCONVEX: 判断球面多边形凸性
  - GCINOUT: 点是否在球面多边形内
  - GCDIST: 球面两点间最短距离
  - GCONARC: 判断点是否在大圆弧上
  - GCCONVERT: 球面角度单位转换
  - GCINTERP: 大圆插值
  - GCQAREA,GCTAREA: 球面四边形和三角形面积
  - GCPNT2GC: 点到大圆距离
  - GCDANGLE,GCAANGLE: 大圆间夹角
- 输入参数:
  - LAT,LON: 纬度经度(度)
  - NPTS: 点数
  - TYPE: 转换单位类型
- 输出参数:
  - 距离、角度、面积等几何量
- 功能特点:
  - 完整的球面几何计算库
  - 支持度、弧度、米、千米、英尺、英里转换
  - 基于精确的球面三角学公式
  - 自动处理数值精度问题
  - 适用于GIS、测地学、球面网格生成

## 126. shsgc_R42.f [📝 src](fortran/shsgc_R42.f)
**功能**: R42菱形截断球谐合成(硬编码版本)

**详细说明**:
- 包含多个子程序:
  - DSHSGCR42: R42球谐合成主接口
  - DTRME: 球谐系数到网格场变换
  - DALPMN1: 关联勒让德函数计算
  - DCOOL,DFIXRL,DGAUSSV等: FFT和高斯网格支持函数
- 输入参数:
  - A,B: 球谐系数实部和虚部(43×43)
  - WORK: 工作数组(LWORK)
- 输出参数:
  - GVAR: 网格变量(128×108)
- 功能特点:
  - 硬编码R42菱形截断参数(IRES=42)
  - 固定128经度×108高斯纬度网格
  - 基于中科院包青提供的算法
  - 包含完整的FFT变换和高斯积分支持
  - 使用预计算的勒让德函数提高效率
  - 适用于T42全球谱模式的网格变换

## 127. simpson_dp.f [📝 src](fortran/simpson_dp.f)
**功能**: Simpson积分规则(双精度，支持缺失值)

**详细说明**:
- 包含三个子程序:
  - DSIMPEQ: 等间距Simpson积分主接口
  - DSIMPNE: 不等间距Simpson积分(支持缺失值)
  - DSIMPNE2: 不等间距Simpson积分核心算法
- 输入参数:
  - X,Y: 自变量和因变量数组(NUM点)
  - DX: 等间距情况的步长
  - YMSG: Y的缺失值代码
- 输出参数:
  - ANS: 积分结果
- 功能特点:
  - 基于3点拉格朗日插值多项式的解析积分
  - 等间距: 误差$\propto$h$^4$f⁽$^4$⁾，h为步长
  - 不等间距: 误差$\propto$最大相邻间距的4次方×f⁽$^4$⁾
  - 自动跳过缺失值对
  - 前端驱动程序简化接口
  - 基于NCAR NSSL算法修改版本

## 128. skewT.f [📝 src](fortran/skewT.f)
**功能**: Skew-T/log-P图热力学计算和绘图工具

**详细说明**:
- 包含多个热力学函数:
  - DSKEWTY: 压力到Skew-T坐标Y轴转换
  - DSKEWTX: 温度到Skew-T坐标X轴转换
  - DTMRSKEWT: 混合比线温度计算
  - DTDASKEWT: 干绝热线温度计算
  - DSATLFTSKEWT: 湿绝热线温度计算
  - DPTLCLSKEWT: 抬升凝结高度计算
  - DSHOWALSKEWT: Showalter稳定指数
  - DPWSKEWT: 可降水量计算
  - DCAPETHERMO: CAPE对流有效位能
- 坐标转换公式:
  - $Y = 132.182 - 44.061 \times \log_{10}(P)$
  - $X = 0.54 \times T + 0.90692 \times Y$
- 功能特点:
  - 完整的大气热力学诊断工具
  - 基于Stipanuk(1973)和标准气象公式
  - 支持干湿绝热过程计算
  - 适用于大气稳定度分析、对流预报

## 129. slatec_DPOLFT.f [📝 src](fortran/slatec_DPOLFT.f)
**功能**: SLATEC多项式最小二乘拟合(双精度)

**详细说明**:
- 包含三个主要子程序:
  - POLFT/POLFTMSG: NCL接口封装程序
  - DPOLFT: SLATEC核心多项式拟合算法
  - DP1VLU: 多项式及导数求值
  - DPCOEF: 正交多项式系数转换为幂级数
- 输入参数:
  - X,Y,W: 数据点、函数值、权重数组(N点)
  - MAXDEG: 最大拟合次数
  - EPS: 拟合精度控制参数
- 输出参数:
  - COEF: 多项式系数(幂级数形式)
  - NDEG: 实际拟合次数
  - R: 拟合值
- 功能特点:
  - 基于正交多项式的稳定算法
  - 支持加权最小二乘拟合
  - 三种次数选择模式: 统计检验、固定次数、误差控制
  - F检验显著性水平: 0.01、0.05、0.10
  - 自动处理缺失值
  - 数值稳定的递推关系

## 130. smth9_dp.f [📝 src](fortran/smth9_dp.f)
**功能**: 9点二维平滑滤波器

**详细说明**:
- 包含一个子程序:
  - DSMTH9: 9点加权平滑核心算法
- 输入参数:
  - X: 输入/输出二维数组(NI×NJ)
  - WRK: 工作数组(NI×NJ)
  - P,Q: 平滑权重参数(建议P=0.5，Q=0.25)
  - XMSG: 缺失值代码
  - LWRAP: 边界环绕标志
- 平滑公式:
  - $f_0' = f_0 + \frac{P}{4}\left(f_2 + f_4 + f_6 + f_8 - 4f_0\right) + \frac{Q}{4}\left(f_1 + f_3 + f_5 + f_7 - 4f_0\right)$
  - 9点模板: 中心点周围8个相邻点
- 功能特点:
  - 基于Olson修改的9点平滑算法
  - 支持周期性边界条件(经度方向)
  - 自动跳过缺失值点
  - 仅平滑内部点，保持边界不变
  - 适用于气象场平滑、噪声滤除

## 131. sp2_util.f [📝 src](fortran/sp2_util.f)
**功能**: SPHEREPACK第2版实用工具库(单精度)

**详细说明**:
- 包含多个球面坐标变换工具:
  - GEOMATV: 地理坐标向量场到数学坐标变换
  - MATGEOV: 数学坐标向量场到地理坐标变换
  - GEOMAT: 地理坐标标量场到数学坐标变换
  - MATGEO: 数学坐标标量场到地理坐标变换
  - TPLATX: 数组转置工具
  - CLATX: 纬度数组翻转工具
  - TRCTPR: 球谐系数截断和锥化
  - GEOSCL: 数组标量乘法
- 坐标变换特征:
  - 地理坐标: 经度×纬度，v分量向北为正
  - 数学坐标: 纬度×经度，v分量取负
  - 包含数组维度变换和纬度翻转
- 功能特点:
  - 基于John Adams原始代码
  - 完整的球面坐标系统转换
  - 支持向量和标量场变换
  - 包含球谐截断和指数锥化功能
  - 适用于SPHEREPACK球谐分析

## 132. sp2_util_dp.f [📝 src](fortran/sp2_util_dp.f)
**功能**: SPHEREPACK第2版实用工具库(双精度版本)

**详细说明**:
- 包含与sp2_util.f相同的功能，全部双精度版本:
  - DGEOMATV: 地理到数学坐标向量变换
  - DMATGEOV: 数学到地理坐标向量变换
  - DGEOMAT: 地理到数学坐标标量变换
  - DMATGEO: 数学到地理坐标标量变换
  - DTPLATX: 双精度数组转置
  - DCLATX: 双精度纬度翻转
  - DTRCTPR: 双精度球谐截断
  - DGEOSCL: 双精度数组缩放
- 变换公式:
  - 地理→数学: 转置+纬度翻转+v取负
  - 数学→地理: v取负+纬度翻转+转置
- 功能特点:
  - 与单精度版本算法完全一致
  - 提供更高的数值精度
  - 适用于高精度球谐分析计算
  - 支持大规模全球模式数据处理

## 133. spareapoly.f [📝 src](fortran/spareapoly.f)
**功能**: 球面多边形面积计算(Bevis-Cambareri算法)

**详细说明**:
- 包含三个子程序:
  - SPAREAPOLYI: 球面多边形面积计算接口
  - SPAREAPOLY: 球面多边形面积核心算法
  - TRNSFRMLON: 球面坐标变换辅助函数
- 输入参数:
  - VLAT,VLON: 多边形顶点纬度、经度(度，NV个)
  - RAD: 球体半径
- 输出参数:
  - AREA: 球面多边形面积(RAD$^2$单位)
- 理论基础:
  - 基于Bevis & Cambareri (1987)算法
  - 面积公式: $A = \left(\sum\text{内角} - \pi(n-2)\right) \times R^2$
  - 顶点按顺时针方向编号
- 功能特点:
  - 适用于任意形状球面多边形
  - 自动处理重复顶点
  - 支持凸多边形和凹多边形
  - 纬度$\pm$90°，经度$\pm$180°
  - 适用于GIS面积计算、地理分析

## 134. specx_dp.f [📝 src](fortran/specx_dp.f)
**功能**: 功率谱密度估计(双精度，完整统计分析)

**详细说明**:
- 包含三个主要子程序:
  - DSPECX: 单变量功率谱分析主接口
  - DSPECXY: 双变量交叉谱分析主接口
  - DSPECXD: 谱分析核心驱动程序
  - DTAPERX: 余弦钟形锥化窗函数
  - DSWRNAV: 谱估计平滑处理
- 输入参数:
  - X,Y: 时间序列数据(NX点)
  - IOPT: 去趋势选项(0=去均值,1=线性,2=二次)
  - JAVE: 谱估计平滑点数(奇数)
  - PCT: 锥化比例[0-1]
  - SCL: 标准化选项(-1,0,1,2或用户定义)
- 输出参数:
  - FRQ,SPCX: 频率和功率谱
  - COSPC,QUSPC: 协谱和正交谱
  - COHER,PHASE: 相干性和相位谱
  - SINFO: 统计信息(自由度、置信区间等)
- 功能特点:
  - 基于FFT的周期图方法
  - 修正Daniell平滑算法
  - 完整的统计显著性检验
  - 支持余弦钟形锥化窗
  - 适用于气候时间序列频域分析

## 135. spi3.f [📝 src](fortran/spi3.f)
**功能**: 标准化降水指数SPI-3计算(NCDC版本)

**详细说明**:
- 包含多个统计分析子程序:
  - SPI3NCDC: NCL接口主程序
  - SPIPE3: SPI核心计算算法
  - SAMLMR: L矩参数估计
  - PELPE3: Pearson Type III参数估计
  - CDFPE3: Pearson Type III分布函数
  - QUASTN: 标准正态分位数函数
- 输入参数:
  - PP: 月降水数据(NTIM，从1月开始)
  - NRUN: 滑动积累月数(典型值3,6,12,24)
  - PMSG: 缺失值代码
- 输出参数:
  - SPI3: 标准化降水指数
  - PROBNE: 累积概率
  - PCPACC: 累积降水量
- 算法特征:
  - 基于L矩的Pearson III型分布拟合
  - 混合分布处理零降水频率
  - SPI范围约束在$\pm$3.09内
  - 按月分别计算分布参数保持季节性
- 功能特点:
  - 美国NCDC官方SPI算法
  - 适用于干旱监测和气候分析
  - 支持多时间尺度(3,6,12,24月等)
  - 包含完整统计检验和L矩估计

## 136. spid.f [📝 src](fortran/spid.f)
**功能**: 标准化降水指数SPI计算(伽马分布版本)

**详细说明**:
- 包含多个统计分析子程序:
  - SPIGAMD: SPI伽马分布主接口
  - ANVNRMD: 标准正态逆变换函数
  - GAMFITD: 伽马分布参数估计
  - GAMCDFD: 伽马分布累积概率函数
  - GAMINVD: 伽马分布逆函数
  - 不完全伽马函数系列: GAMMAPD, GAMMAQD, GAMSERD, GAMMCFD等
- 输入参数:
  - PP: 降水时间序列数据(NTIM)
  - NRUN: 累积时间尺度(月数)
  - PMSG: 缺失值代码
- 输出参数:
  - INDEX: 标准化降水指数值
- 数学原理:
  - 伽马分布概率密度: $f(x) = \frac{1}{\beta^\gamma \Gamma(\gamma)} x^{\gamma-1} e^{-x/\beta}$
  - 参数估计: $\alpha = \ln(\bar{x}) - \frac{\sum \ln(x)}{n}$, $\gamma = \frac{1 + \sqrt{1 + 4\alpha/3}}{4\alpha}$
  - 标准化变换: $\text{SPI} = \Phi^{-1}[F_\gamma(x)]$
- 功能特点:
  - 基于Thom(1958)极大似然估计方法
  - 支持含零值的混合伽马分布
  - 使用Abramowitz & Stegun正态逆变换
  - 按月分别拟合分布参数
  - 适用于干旱监测和降水异常分析

## 137. stattrop_dp.f [📝 src](fortran/stattrop_dp.f)
**功能**: WMO对流层高度计算(双精度)

**详细说明**:
- 包含两个子程序:
  - STATTROPX: NCL接口封装程序
  - DSTATTROP: WMO对流层高度核心算法
- 输入参数:
  - PFULL,TFULL: 层压力和温度(NLEV层)
  - LAPSEC: 温度递减率阈值(典型值2 K/km)
  - PUNIT: 压力单位(0=hPa，1=Pa)
- 输出参数:
  - PTROP: 对流层顶压力
  - ITROP: 对流层顶层索引
- WMO定义准则:
  - 第一对流层顶: 温度递减率$\leq$2 K/km的最低层
  - 2公里平均检验: 该层上方2km内平均递减率$\leq$2 K/km
- 计算方法:
  - 温度递减率: $\gamma = -\frac{dT}{dz} = \frac{d\ln T}{d\ln p} \times \frac{g}{R} \times 1000$
  - 静力学近似: $dz = -\frac{dp}{g\rho} = -\frac{RT}{gp}dp$
  - 对数插值: $p_{trop} = \exp[p_1 + w(p_2 - p_1)]$
- 功能特点:
  - 基于WMO(1992)国际气象词汇定义
  - 支持450-85 hPa范围对流层顶识别
  - 线性对数压力插值提高精度
  - 兼容NCEP再分析对流层顶算法
  - 适用于大气垂直结构分析

## 138. statx.f [📝 src](fortran/statx.f)
**功能**: 统计分析工具库(单精度，支持缺失值)

**详细说明**:
- 包含完整的统计分析子程序集:
  - STAT2/STAT4: 二阶/四阶矩估计
  - STAT2T: 截尾均值和方差估计
  - MEDMRNG: 中位数、极差、中程估计
  - ESAUTO/ESCROS: 自相关和交叉相关函数
  - VCVMNS/CORMNS: 方差协方差和相关矩阵
  - ERRCHK: 异常值检测
  - SORTU/ISORTU: 快速排序算法
- 统计量计算:
  - 样本方差: $s^2 = \frac{\sum(x_i - \bar{x})^2}{n-1}$
  - 偏度系数: $\text{Skew} = \frac{\sum(x_i - \bar{x})^3}{ns^3}$
  - 峰度系数: $\text{Kurt} = \frac{\sum(x_i - \bar{x})^4}{ns^4} - 3$
  - 自相关函数: $r(\tau) = \frac{\sum(x_t - \bar{x})(x_{t+\tau} - \bar{x})}{\sum(x_t - \bar{x})^2}$
- 功能特点:
  - 全面的缺失值处理机制
  - 截尾统计降低极值影响
  - 数值稳定的递推算法
  - 对称存储模式节省内存
  - 包含完整的矩阵统计运算
  - 适用于时间序列和空间数据分析

## 139. statx_dp.f [📝 src](fortran/statx_dp.f)
**功能**: 统计分析工具库(双精度版本)

**详细说明**:
- 包含与statx.f相同功能的双精度实现:
  - DSTAT2/DSTAT4: 双精度矩估计
  - DSTAT2T: 双精度截尾统计
  - DMEDMRNG: 双精度中位数计算
  - DESAUTO/DESCROS: 双精度相关分析
  - DVCVMNS/DCORMNS: 双精度矩阵统计
  - DSORTU/DISORTU: 双精度排序算法
  - COLLAPSEXY: 缺失值配对数据提取
- 数值精度优势:
  - 双精度浮点运算(约15位有效数字)
  - 减少累积舍入误差
  - 提高大样本统计的数值稳定性
  - 支持高精度科学计算需求
- 功能特点:
  - 与单精度版本算法完全一致
  - 保持相同的接口规范
  - 适用于高精度统计分析
  - 大数据集处理的数值可靠性
  - 科学计算和工程应用的首选

## 140. stndAtmUS76_dp.f [📝 src](fortran/stndAtmUS76_dp.f)
**功能**: 美国标准大气1976版热力学计算

**详细说明**:
- 包含三个子程序:
  - DSTDATMZ: 高度→温度、密度、压力
  - DSTDATMP: 压力→高度、温度、密度  
  - STDZ2P/STDP2Z: 标准大气核心算法
- 标准大气层结构(0-84.852km):
  - 对流层(0-11km): $T = 15.0 - 6.5z$ (°C)
  - 平流层下层(11-20km): $T = -56.5$ (°C)  
  - 平流层中层(20-32km): $T = -56.5 + 1.0(z-20)$ (°C)
  - 平流层上层(32-47km): $T = -44.5 + 2.8(z-32)$ (°C)
- 热力学关系:
  - 静力学方程: $dp = -\rho g dz$
  - 理想气体: $p = \rho R T$, $R = 287.04$ J/(kg·K)
  - 等温层: $p = p_0 \exp\left(-\frac{gz}{RT}\right)$
  - 线性梯度层: $p = p_0 \left(\frac{T}{T_0}\right)^{-\frac{g}{R\gamma}}$
- 功能特点:
  - 基于美国标准大气1976版
  - 有效高度范围0-84.852km
  - 压力范围1013.25-0.00373hPa
  - 支持米和英尺单位转换
  - 高精度双精度计算
  - 适用于航空航天和大气物理应用

## 141. student_tprob.f [📝 src](fortran/student_tprob.f)
**功能**: 学生t分布概率计算(基于不完全贝塔函数)

**详细说明**:
- 包含一个子程序:
  - STUPROBT: 学生t分布双尾概率计算
- 输入参数:
  - T: t统计量值
  - DF: 自由度(>0)
  - TMSG: 缺失值代码
- 输出参数:
  - RESULT: 双尾概率值P(|T|>t)
- 数学原理:
  - t分布概率密度: $f(t) = \frac{\Gamma(\frac{\nu+1}{2})}{\sqrt{\nu\pi}\Gamma(\frac{\nu}{2})} \left(1+\frac{t^2}{\nu}\right)^{-\frac{\nu+1}{2}}$
  - 累积概率变换: $P(|T|>t) = I_{\frac{\nu}{\nu+t^2}}\left(\frac{\nu}{2}, \frac{1}{2}\right)$
  - 其中$I_x(a,b)$为不完全贝塔函数
- 功能特点:
  - 基于SLATEC库不完全贝塔函数DBETAI
  - 通过变换避免直接计算t分布积分
  - 自动处理缺失值和边界条件
  - 适用于假设检验和置信区间计算

## 142. svd_lap_dp.f [📝 src](fortran/svd_lap_dp.f)
**功能**: LAPACK奇异值分解(双精度，基于Bretherton方法)

**详细说明**:
- 包含三个主要子程序:
  - DSVDLAP: SVD主接口(完整统计分析)
  - DSVDSV: SVD简化接口(仅返回奇异向量)  
  - DSVDLAPACK1/DSVDLAPSV: SVD核心算法
- 输入参数:
  - X,Y: 输入数据矩阵(LDXY×NCX, LDXY×NCY)
  - MRT: 时间维度长度
  - NCX,NCY: X和Y的空间维度
  - NSV: 需要的奇异值个数
  - IFLAG: 标准化选项(1=标准化)
- 输出参数:
  - HOMLFT/HOMRGT: 左/右同类相关图
  - HETLFT/HETRGT: 左/右异类相关图
  - AK,BK: 时间展开系数
  - PCVAR: 方差贡献百分比
- SVD分解过程:
  - 交叉协方差矩阵: $C = \frac{1}{T}X^T Y$
  - SVD分解: $C = U \Sigma V^T$
  - 时间展开系数: $A_k(t) = X(t) \cdot U_k$, $B_k(t) = Y(t) \cdot V_k$
- 功能特点:
  - 基于LAPACK DGESVD高精度SVD算法
  - 完整的Bretherton等人统计诊断
  - 支持数据标准化和缺失值处理
  - 包含条件数和秩估计
  - 适用于耦合模态分析和降维

## 143. svd_util_dp.f [📝 src](fortran/svd_util_dp.f)
**功能**: SVD分析实用工具集(双精度)

**详细说明**:
- 包含四个实用子程序:
  - DSVDMTR: 矩阵转置工具
  - DSVDINFO: SVD信息统计分析
  - DSVDPAR: 矩阵数据打印输出
  - DSVDPVC: 向量数据打印输出
- DSVDINFO功能:
  - 计算Frobenius范数: $||A||_F = \sqrt{\sum_{i=1}^{\min(m,n)} \sigma_i^2}$
  - 条件数估计: $\kappa = \frac{\sigma_{\max}}{\sigma_{\min}}$
  - 数值秩确定: $\text{rank}(A) = \#\{\sigma_i : \sigma_i > \text{TOL}\}$
  - 方差贡献: $\text{PC}_i = \frac{\sigma_i^2}{\sum \sigma_j^2} \times 100\%$
- 数值稳定性:
  - 基于机器精度的容差设定
  - 自动处理奇异值排序
  - 稳健的秩估计算法
- 功能特点:
  - 为SVD分析提供完整的辅助工具
  - 标准化的输出格式便于调试
  - 双精度数值精度保证
  - 模块化设计便于维护

## 144. taper.f [📝 src](fortran/taper.f)
**功能**: 分段余弦钟形锥化窗函数

**详细说明**:
- 包含一个子程序:
  - DTAPER: 时间序列锥化处理
- 输入参数:
  - X: 输入时间序列(N点)
  - P: 锥化比例[0,1](如P=0.1表示10%)
  - IOPT: 锥化选项(0=锥化到序列均值，1=强制锥化到0)
- 输出参数:
  - XT: 锥化后的时间序列
- 锥化窗函数:
  - 锥化点数: $M = \max(1, \lfloor P \times N/2 + 0.5 \rfloor)$
  - 权重函数: $w(i) = 0.5 - 0.5\cos\left(\frac{\pi(i-0.5)}{M}\right)$, $i = 1,2,\ldots,M$
  - 锥化公式: $x_t'(i) = (x(i) - \bar{x}) \times w(i) + \bar{x}$
- 应用位置:
  - 序列头部: $i = 1, 2, \ldots, M$
  - 序列尾部: $i = N-M+1, N-M+2, \ldots, N$
  - 中间部分保持不变
- 功能特点:
  - 基于Bloomfield傅里叶分析理论
  - 分段余弦钟形平滑过渡
  - 自动检测常数序列并强制锥化到0
  - 减少FFT分析中的频谱泄漏
  - 适用于非周期信号的频域分析

## 145. tdez1d.f [📝 src](fortran/tdez1d.f)
**功能**: NCAR Graphics三维轨迹标记绘图

**详细说明**:
- 包含一个绘图子程序:
  - TDEZ1D: 三维空间轨迹标记可视化
- 输入参数:
  - X,Y,Z: 三维轨迹坐标数组(NX点)
  - IMRK: 标记类型(1-5: 四面体、八面体、立方体、二十面体、球体)
  - RMRK: 标记半径
  - SMRK: 标记间隔(负值允许重叠)
  - RMULT,THETA,PHI: 观察点球坐标(倍数、方位角、仰角)
  - IST: 着色风格索引(1-8种颜色方案)
- 观察点计算:
  - 默认位置: $(R, \theta, \phi) = (2.5 \times DL, -55°, 70°)$
  - 笛卡尔坐标转换: 
    - $x_{eye} = x_{mid} + R \sin\phi \cos\theta$
    - $y_{eye} = y_{mid} + R \sin\phi \sin\theta$ 
    - $z_{eye} = z_{mid} + R \cos\phi$
  - 其中$DL = \sqrt{(\Delta x)^2 + (\Delta y)^2 + (\Delta z)^2}$为对角线长度
- 颜色方案:
  - IST=1: 线框模式
  - IST=2-8: 底部灰度+顶部彩色(红、绿、蓝、青、紫、黄)
  - 负值IST: 白色背景+黑色前景
- 功能特点:
  - 基于NCAR Graphics Tdpack三维绘图库
  - 支持多种几何标记和着色方案
  - 自动计算最佳观察角度和标记尺寸
  - 三角剖分渲染实现真实感显示
  - 适用于三维轨迹和散点数据可视化

## 146. thornthwaite_v2.f [📝 src](fortran/thornthwaite_v2.f)
**功能**: Thornthwaite蒸散发潜力计算(第2版，支持缺失值)

**详细说明**:
- 包含一个子程序:
  - THORN2: Thornthwaite蒸散发潜力计算主接口
- 输入参数:
  - T: 月平均温度数组(NTIM，°C)
  - LAT: 纬度(度)
  - TMSG: 缺失值代码
- 输出参数:
  - ETP: 蒸散发潜力(NTIM，mm/月)
- Thornthwaite公式:
  - 热量指数: $J = \sum_{i=1}^{12} \left(\frac{T_i}{5}\right)^{1.514}$ (仅当$T_i > 0$)
  - 指数系数: $c = 6.75 \times 10^{-7} J^3 - 7.71 \times 10^{-5} J^2 + 1.792 \times 10^{-2} J + 0.49239$
  - 蒸散发: $\text{ETP} = K \times 16 \times \left(\frac{10T}{J}\right)^c$ (当$T > 0$时)
- 日照修正系数:
  - 太阳时角: $\omega = \arccos(-\tan\phi \tan\delta)$
  - 日照时数: $N = \frac{24\omega}{\pi}$
  - 修正系数: $K = \frac{N}{12} \times \frac{\text{days}}{30}$
- 功能特点:
  - 基于原始C代码http://sac.csic.es/spei/spei_index.html
  - 使用月平均太阳赤纬角预计算表
  - 自动处理缺失值和负温度
  - 考虑纬度和季节变化的日照修正
  - 适用于气候学蒸散发估算和干旱指数计算

## 147. triple2grid2d.f [📝 src](fortran/triple2grid2d.f)
**功能**: 散点数据到规则网格的最近邻插值

**详细说明**:
- 包含一个子程序:
  - TRIPLE2GRID2D: 三元组数据到二维网格插值
- 输入参数:
  - X,Y,Z: 观测点坐标和数值(KPTS个观测点)
  - LAT2D,LON2D: 目标网格纬度和经度(MDIM×NDIM)
  - DISTMX: 最大搜索距离
  - MOPT: 距离计算选项(0=直角坐标，1=大圆距离)
- 输出参数:
  - ZGRID: 插值后的网格数据(MDIM×NDIM)
- 距离计算公式:
  - 直角坐标: $d = \sqrt{(x_i - x_g)^2 + (y_i - y_g)^2}$
  - 大圆距离: $d = R_e \arccos[\sin\phi_g \sin\phi_i + \cos\phi_g \cos\phi_i \cos(\lambda_i - \lambda_g)]$
  - 其中$R_e = 6371.22$ km为地球半径
- 功能特点:
  - 最近邻插值算法：每个网格点赋值为距离最近的观测值
  - 支持直角坐标和球面坐标两种距离计算
  - 距离阈值控制：超出DISTMX的观测点不参与插值
  - 自动处理缺失值传播
  - 适用于气象站数据到网格的快速插值

## 148. trssphx.f [📝 src](fortran/trssphx.f)
**功能**: 球谐函数标量场网格转换(SPHEREPACK接口)

**详细说明**:
- 包含一个子程序:
  - TRSSPHX: 标量场球面网格间转换主接口
- 输入参数:
  - INTL: 初始化标志(0=初始化，1=转换)
  - IGRIDA/IGRIDB: 源/目标网格类型([$\pm$1,0/1]：等间距/高斯，纬度/经度优先)
  - NLONA,NLATA/NLONB,NLATB: 源/目标网格维度
  - DA: 源网格标量场数据
  - JMODE: 截断和锥化选项
- 输出参数:
  - DB: 目标网格标量场数据
- 球谐变换过程:
  - 分析: $f(\lambda,\theta) \rightarrow \{A_n^m, B_n^m\}$
  - 截断: 球谐系数截断到指定波数
  - 锥化: 可选的指数锥化减少Gibbs现象
  - 合成: $\{A_n^m, B_n^m\} \rightarrow f'(\lambda',\theta')$
- 功能特点:
  - 基于NCAR SPHEREPACK双精度库
  - 支持等间距和高斯网格互转
  - 支持不同分辨率网格变换
  - 可选三角截断和指数锥化
  - 保持工作空间优化的高效算法
  - 适用于全球模式数据后处理

## 149. trvsphx.f [📝 src](fortran/trvsphx.f)
**功能**: 球谐函数向量场网格转换(SPHEREPACK接口)

**详细说明**:
- 包含一个子程序:
  - TRVSPHX: 向量场球面网格间转换主接口
- 输入参数:
  - INTL: 初始化标志(0=初始化，1=转换)
  - IGRIDA/IGRIDB: 源/目标网格类型
  - UA,VA: 源网格向量场分量
  - IVECA/IVECB: 向量场类型(0=数学坐标，1=地理坐标)
  - JMODE: 截断和锥化选项
- 输出参数:
  - UB,VB: 目标网格向量场分量
- 向量球谐变换:
  - 分析: $(u,v) \rightarrow \{A_n^m, B_n^m, C_n^m, D_n^m\}$
  - 其中$(A,B)$为涡度位势，$(C,D)$为散度位势
  - 截断和锥化处理应用于所有系数
  - 合成: $\{A_n^m, B_n^m, C_n^m, D_n^m\} \rightarrow (u',v')$
- 坐标系统转换:
  - 地理坐标：u向东，v向北
  - 数学坐标：u向东，v向南(负向)
  - 自动处理纬度排序(北→南或南→北)
- 功能特点:
  - 基于NCAR SPHEREPACK向量球谐库
  - 保持向量场散度和涡度性质
  - 支持坐标系统和网格类型转换
  - 可选截断和锥化优化
  - 适用于风场、海流等向量场数据处理

## 150. ttest.f [📝 src](fortran/ttest.f)
**功能**: 统计假设检验工具集(t检验、F检验、相关检验)

**详细说明**:
- 包含四个统计检验子程序:
  - DTTEST: 双样本t检验(等/非等方差)
  - DFTEST: 双样本F检验(方差相等性)
  - DRTEST: 相关系数显著性检验
  - DEQVSIZ: 等效样本量估计(红噪声过程)
- t检验公式:
  - 等方差: $t = \frac{\bar{x}_1 - \bar{x}_2}{\sqrt{s_p^2(1/n_1 + 1/n_2)}}$, $s_p^2 = \frac{(n_1-1)s_1^2 + (n_2-1)s_2^2}{n_1+n_2-2}$
  - 非等方差(Welch): $t = \frac{\bar{x}_1 - \bar{x}_2}{\sqrt{s_1^2/n_1 + s_2^2/n_2}}$, $\nu = \frac{(s_1^2/n_1 + s_2^2/n_2)^2}{(s_1^2/n_1)^2/(n_1-1) + (s_2^2/n_2)^2/(n_2-1)}$
- F检验公式:
  - $F = \frac{s_{\max}^2}{s_{\min}^2}$，自由度$\nu_1 = n_{\max} - 1$, $\nu_2 = n_{\min} - 1$
- 相关检验:
  - $t = r\sqrt{\frac{n-2}{1-r^2}}$，自由度$\nu = n-2$
  - 等效样本量: $n_{eff} = n \frac{1-r_1}{1+r_1}$ (当红噪声显著时)
- 功能特点:
  - 基于SLATEC库不完全贝塔函数计算p值
  - 完整的假设检验统计量和显著性水平
  - 支持红噪声时间序列的等效自由度修正
  - 涵盖气候统计中最常用的检验方法
  - 适用于气象数据的统计显著性分析

## 151. upsl_dp.f [📝 src](fortran/upsl_dp.f)
**功能**: 使用水位测高方程计算海平面压力

**详细说明**:
- 包含两个子程序:
  - DPSLHY1: 使用虚温的海平面压力计算
  - DPSLHY2: 使用温度和混合比的海平面压力计算
- 输入参数:
  - PRES: 最低垂直层压力(Pa)
  - Z: 最低垂直层位势高度(m)
  - TV: 虚温(K) 或 T: 温度(K)
  - W: 混合比(kg/kg)(DPSLHY2)
  - XMSG: 缺失值代码
- 输出参数:
  - DPSLHY1/DPSLHY2: 海平面压力(Pa)
- 功能特点:
  - 基于水位测高方程: $p_{sl} = p \cdot \exp\left(\frac{gz}{R \cdot T_v}\right)$
  - DPSLHY2中虚温计算: $T_v = T \cdot (1 + 0.608w)$
  - 使用标准重力加速度: $g = 9.80665$ m/s$^2$
  - 干空气气体常数: $R = 287.04$ J/(kg·K)
  - 自动处理缺失值
  - 适用于气象站观测数据的海平面压力订正

## 152. varimax_JiangLing_dp.f [📝 src](fortran/varimax_JiangLing_dp.f)
**功能**: 改进的Varimax旋转算法,用于EOF分析的因子旋转

**详细说明**:
- 包含两个子程序:
  - ZROTEOF: 主旋转接口程序,处理缺失值检测和方差计算
  - ZVORS: 核心Varimax旋转算法(改进的收敛判据)
- 输入参数:
  - V: 特征向量矩阵(ND×NF)
  - NV,NF,ND: 变量数、因子数、数组维度
  - EVAL: 特征值数组(NF)
  - VMSG: 缺失值代码
  - IOPT: 缩放选项(0=不缩放,1=缩放)
- 输出参数:
  - V: 旋转后的特征向量矩阵
  - ROTPCV: 旋转后方差贡献百分比(NF)
  - ROTVAR: 旋转后方差(NF)
  - KFLAG: 缺失值存在标志
- 功能特点:
  - 基于Kaiser行标准化的Varimax准则
  - 改进的收敛判据: $|\gamma(k+1) - \gamma(k)| < 10^{-5}$
  - 旋转准则: $\gamma = \sum_{j=1}^p \left[\sum_{i=1}^n a_{ij}^4 - \frac{1}{n}\left(\sum_{i=1}^n a_{ij}^2\right)^2\right]$
  - 支持缺失值的稳健处理
  - 适用于气候模态分析的因子简单结构识别

## 153. varimax_dp.f [📝 src](fortran/varimax_dp.f)
**功能**: 经典Varimax因子旋转算法,用于主成分分析

**详细说明**:
- 包含四个子程序:
  - ROTEOF: 主旋转接口程序
  - VORS: 核心Varimax旋转算法
  - DSUMFR: 向量求和/平方和计算
  - DSCPFR: 向量标量积计算
  - VORSMSG: 缺失值处理版本旋转
- 输入参数:
  - V: 待旋转矩阵(ND×NF)
  - EVAL: 特征值(NF)
  - IOPT: 处理选项($\leq$0=不缩放,>0=缩放)
- 输出参数:
  - V: 旋转后的因子载荷矩阵
  - ROTPCV: 旋转后方差贡献百分比
  - ROTVAR: 旋转后方差
- 功能特点:
  - 基于Veldman(1967)经典Varimax算法
  - Kaiser行标准化: $v_{ij}' = \frac{v_{ij}}{\sqrt{\sum_k v_{ik}^2}}$
  - 旋转角度: $\tan(4\theta) = \frac{D - 2AB/n}{C - (A^2-B^2)/n}$
  - 迭代收敛精度: $\epsilon = 0.0001$ 弧度
  - 兼容IMSL FROTA函数(w=1.0, norm=1, maxit=30)
  - 适用于因子分析简单结构获取

## 154. vibeta_dp.f [📝 src](fortran/vibeta_dp.f)
**功能**: 使用Kuo/Trenberth诊断方法进行垂直积分

**详细说明**:
- 包含四个子程序:
  - DVIBETA: 垂直积分主算法
  - DINT2P2: 压力坐标插值
  - DINITVC: 向量初始化
  - DCOPYVC: 向量复制
- 输入参数:
  - P: 压力层数组(NLEV,Pa)
  - X: 待积分量数组(NLEV)
  - PSFC: 地面压力(Pa)
  - PBOT,PTOP: 积分下界和上界压力(Pa)
  - LINLOG: 插值类型(1=线性,2=对数)
  - XMSG: 缺失值代码
- 输出参数:
  - VINT: 垂直积分结果
  - IER: 错误代码
- 功能特点:
  - 基于Kuo & Trenberth(1991)诊断方法
  - 积分公式: $\int_{P_{top}}^{P_{bot}} X \frac{dp}{g} = \sum_{k} \beta_k X_k \Delta p_k$
  - 权重系数: $\beta_k = \frac{p_{sfc} - p_{k+1}}{p_{k-1} - p_{k+1}}$ (表面层)
  - 支持缺失值的探空数据处理
  - 适用于大气柱总量计算(水汽、臭氧等)

## 155. vinth2p_dp.f [📝 src](fortran/vinth2p_dp.f)
**功能**: CCM2/3混合坐标数据到压力坐标的垂直插值

**详细说明**:
- 包含两个子程序:
  - VINTH2P: 主垂直插值程序
  - PRNT: 调试打印工具
- 输入参数:
  - DATI: 混合坐标数据(IMAX×NLAT×NLEVI)
  - HBCOFA,HBCOFB: 混合坐标系数A和B(NLEVIP1)
  - P0: 基准压力(mb)
  - PLEVO: 输出压力层(NLEVO,mb)
  - PSFC: 地面压力(IMAX×NLAT,Pa)
  - INTYP: 插值类型(1=线性,2=对数,3=双对数)
  - KXTRP: 外推标志(0=用特殊值,1=外推)
- 输出参数:
  - DATO: 压力坐标数据(IMAX×NLAT×NLEVO)
- 功能特点:
  - 混合坐标公式: $P(k) = A(k) \times P_0 + B(k) \times P_{sfc}$
  - 线性插值: $f = f_1 + (f_2-f_1) \frac{p-p_1}{p_2-p_1}$
  - 对数插值: $f = f_1 + (f_2-f_1) \frac{\ln(p/p_1)}{\ln(p_2/p_1)}$
  - 双对数插值: $f = f_1 + (f_2-f_1) \frac{A2LN(p)-A2LN(p_1)}{A2LN(p_2)-A2LN(p_1)}$
  - 其中$A2LN(x) = \ln(\ln(x+2.72))$
  - 支持边界外推和特殊值处理
  - 适用于大气模式垂直坐标转换

## 156. vinth2p_ecmwf.f [📝 src](fortran/vinth2p_ecmwf.f)
**功能**: ECMWF模式混合坐标数据到压力坐标的垂直插值(含ECMWF外推)

**详细说明**:
- 包含一个子程序:
  - VINTH2PECMWF: 主垂直插值程序(带ECMWF外推算法)
- 输入参数:
  - DATI: 混合坐标数据(IMAX×NLAT×NLEVI)
  - HBCOFA,HBCOFB: 混合坐标系数A和B(NLEVIP1)
  - P0: 基准压力(mb)
  - PLEVO: 输出压力层(NLEVO,mb)
  - PSFC: 地面压力(IMAX×NLAT,Pa)
  - VARFLG: 变量标志(-1=位势高度,+1=温度,0=其他)
  - TBOT: 最低层温度(IMAX×NLAT,K)
  - PHIS: 地面位势(IMAX×NLAT,m$^2$/s$^2$)
  - INTYP: 插值类型(1=线性,2=对数,3=双对数)
- 输出参数:
  - DATO: 压力坐标数据(IMAX×NLAT×NLEVO)
- 功能特点:
  - 基于ECMWF IFS模式的外推算法
  - 温度外推: $T = T^* \left(1 + \alpha \ln p + \frac{1}{2}\alpha^2(\ln p)^2 + \frac{1}{6}\alpha^3(\ln p)^3\right)$
  - 位势高度外推: $Z = Z_s - \frac{RT^*}{g}\ln\left(\frac{p}{p_s}\right)\left(1 + \frac{1}{2}\alpha\ln p + \frac{1}{6}\alpha^2(\ln p)^2\right)$
  - 其中$\alpha = \frac{0.0065 \times R}{g}$为标准大气递减率参数
  - 支持复杂的高度修正算法(2000-2500m过渡层)
  - 适用于ECMWF数据的高精度垂直插值

## 157. vinth2p_ecmwf_nodes_dp.f [📝 src](fortran/vinth2p_ecmwf_nodes_dp.f)
**功能**: ECMWF混合坐标节点数据的垂直插值(一维节点数组版本)

**详细说明**:
- 包含一个子程序:
  - DVINTH2PECMWFNODES: 节点版ECMWF垂直插值算法
- 输入参数:
  - DATI: 混合坐标数据(NPTS×NLEVI)
  - HBCOFA,HBCOFB: 混合坐标系数A和B
  - PSFC: 地面压力(NPTS,Pa)
  - TBOT: 最低层温度(NPTS,K)
  - PHIS: 地面位势(NPTS,m$^2$/s$^2$)
  - VARFLG: 变量类型标志
- 输出参数:
  - DATO: 插值后数据(NPTS×NLEVO)
- 功能特点:
  - 适用于一维节点数组的高效处理
  - 使用相同的ECMWF外推算法
  - 温度修正公式: $T^* = T_{bot} \times \left(1 + \alpha\left(\frac{p_{sfc}}{p_{lev}} - 1\right)\right)$
  - 高度相关的温度梯度处理: 2000m以下使用标准梯度,2000-2500m使用渐变
  - 双精度计算保证数值精度
  - 适用于观测点或不规则网格数据插值

## 158. vinth2p_nodes_dp.f [📝 src](fortran/vinth2p_nodes_dp.f)
**功能**: 标准混合坐标节点数据的垂直插值(无ECMWF外推)

**详细说明**:
- 包含一个子程序:
  - DVINTH2PNODES: 标准节点版垂直插值算法
- 输入参数:
  - DATI: 混合坐标数据(NPTS×NLEVI)
  - HBCOFA,HBCOFB: 混合坐标系数
  - PLEVO: 输出压力层(NLEVO,mb)
  - PSFC: 地面压力(NPTS,Pa)
  - KXTRP: 外推标志(0=用特殊值,1=简单外推)
- 输出参数:
  - DATO: 插值后数据(NPTS×NLEVO)
- 功能特点:
  - 不含ECMWF复杂外推算法
  - 三种插值方法: 线性、对数、双对数
  - 双对数函数: $A2LN(x) = \ln(\ln(x + 2.72))$
  - 简单外推: 使用最低层数据值
  - 混合坐标公式: $P(k) = A(k) \times P_0 + B(k) \times P_{sfc}$
  - 适用于标准大气模式数据插值

## 159. vintp2p_ecmwf.f [📝 src](fortran/vintp2p_ecmwf.f)
**功能**: ECMWF模式压力坐标数据到压力坐标的垂直插值

**详细说明**:
- 包含一个子程序:
  - VINTP2PECMWF: 压力到压力插值(含ECMWF外推)
- 输入参数:
  - DATI: 输入压力坐标数据(IMAX×NLAT×NLEVI)
  - PRESI: 输入压力场(IMAX×NLAT×NLEVI,Pa)
  - PLEVO: 输出压力层(NLEVO,mb)
  - PSFC: 地面压力(IMAX×NLAT,Pa)
  - VARFLG: 变量类型(-1=高度,+1=温度,0=其他)
  - TBOT,PHIS: 地面温度和位势
- 输出参数:
  - DATO: 插值后压力坐标数据(IMAX×NLAT×NLEVO)
- 功能特点:
  - 基于已有压力场的插值算法
  - 使用与vinth2p_ecmwf相同的ECMWF外推公式
  - 三种插值类型: $f = f_1 + (f_2-f_1) \times \frac{I(p)-I(p_1)}{I(p_2)-I(p_1)}$
  - 其中$I(p)$为插值函数: 线性、$\ln p$或$\ln(\ln(p+2.72))$
  - 适用于再分析数据和观测数据的垂直插值
  - 不需要混合坐标系数,直接使用压力场

## 160. volAve.f [📝 src](fortran/volAve.f)
**功能**: 计算三维数据场的加权体积平均值

**详细说明**:
- 包含两个子程序:
  - DWGTVOLAVE: 使用分离权重的体积平均
  - DWGTVOLAVECCM: 使用3D权重的体积平均(适用于CCM模式)
- 输入参数:
  - T: 3D数据场(MX×NY×KZ)
  - WGTX: X方向权重(MX)
  - WGTY: Y方向权重(NY)  
  - WGTZ: Z方向权重(KZ)或3D权重(MX×NY×KZ)
  - IFLAG: 缺失值处理标志(0=忽略,1=遇到则返回缺失值)
- 输出参数:
  - AVE: 体积加权平均值
- 功能特点:
  - 两步加权平均算法:
    1. 垂直加权: $\bar{T}_z(x,y) = \frac{\sum_k T(x,y,k) \cdot w_z(k)}{\sum_k w_z(k)}$
    2. 水平加权: $\langle T \rangle = \frac{\sum_{x,y} \bar{T}_z(x,y) \cdot w_x(x) \cdot w_y(y)}{\sum_{x,y} w_x(x) \cdot w_y(y)}$
  - CCM版本支持3D权重(如$\Delta p$场): $w_z(x,y,k) = \Delta p(x,y,k)$
  - 自动处理缺失值
  - 适用于大气柱质量加权平均、全球平均等计算  
## 161. volRmse.f [📝 src](fortran/volRmse.f)
**功能**: 计算两个三维数据场之间的加权体积均方根误差

**详细说明**:
- 包含两个子程序:
  - DWGTVOLRMSE: 使用分离权重的三维RMSE计算
  - DWGTVOLRMSECCM: 使用3D权重的CCM版本RMSE计算
- 输入参数:
  - T,Q: 两个3D数据场(MX×NY×KZ)
  - WGTX: X方向权重(MX)
  - WGTY: Y方向权重(NY) 
  - WGTZ: Z方向权重(KZ)或3D权重(MX×NY×KZ)
  - DPT,DPQ: CCM版本中的3D权重场(如层厚Δp)
  - IFLAG: 缺失值处理标志
- 输出参数:
  - RMSE: 体积加权均方根误差
- 功能特点:
  - 两步RMSE计算:
    1. 垂直RMSE: $\text{RMSE}_z(x,y) = \sqrt{\frac{\sum_k w_z(k) \cdot (T(x,y,k)-Q(x,y,k))^2}{\sum_k w_z(k)}}$
    2. 水平加权: $\text{RMSE} = \sqrt{\frac{\sum_{x,y} w_x(x) \cdot w_y(y) \cdot \text{RMSE}_z^2(x,y)}{\sum_{x,y} w_x(x) \cdot w_y(y)}}$
  - CCM版本使用平均层厚: $w_{avg}(x,y,k) = \frac{1}{2}(w_T(x,y,k) + w_Q(x,y,k))$
  - 适用于三维模式验证、数据对比、模式间误差分析

## 162. vpower.f [📝 src](fortran/vpower.f)
**功能**: 计算两个变量的时空功率谱、交叉谱、相干性和相位

**详细说明**:
- 包含四个子程序:
  - SPCTIMCROSS1: 主控程序，计算数组维度和调用核心算法
  - SPCTIMCROSS2: 核心时空谱分析算法
  - SPCTIMCROSS3: 平滑处理和相干性计算
  - SMTH121STCP: 1-2-1平滑滤波器
- 输入参数:
  - X,Y: 两个3D输入场(NL×NM×NT)，纬度×经度×时间
  - NL,NM,NT: 纬度点数、经度点数、时间点数
- 输出参数:
  - STC: 时空谱结果(NLP1×NT2P1×16)，包含16个谱分量
- 算法流程:
  - 空间FFT: $\tilde{X}(n,t) = \text{FFT}_{\text{lon}}[X(\text{lon},\text{lat},t)]$
  - 时间FFT: $\hat{X}(n,\omega) = \text{FFT}_{\text{time}}[\tilde{X}(n,t)]$
  - 功率谱: $P_{XX}(n,\omega) = |\hat{X}(n,\omega)|^2$
  - 交叉谱: $P_{XY}(n,\omega) = \text{Re}[\hat{X}^*(n,\omega)\hat{Y}(n,\omega)]$
  - 相干性: $\text{Coh}^2(n,\omega) = \frac{P_{XY}^2 + Q_{XY}^2}{P_{XX} \cdot P_{YY}}$
- 功能特点:
  - 区分对称和反对称分量(相对于赤道)
  - 1-2-1平滑滤波减少噪声
  - 计算相位关系和向量分量
  - 适用于大气波动分析、Wheeler-Kiladis图制作

## 163. wavelet.f [📝 src](fortran/wavelet.f)
**功能**: 连续小波变换分析,支持多种母小波和显著性检验

**详细说明**:
- 包含多个子程序:
  - WAVELETI: NCL接口主程序
  - WAVELET: 核心小波变换算法
  - WAVE_FUNCTION: 母小波函数计算
  - WAVE_SIGNIF: 显著性检验
  - STATSWAVE: 时间序列统计量计算
- 支持的母小波:
  - Morlet小波(mother=0): $\psi(t) = \pi^{-1/4}e^{ik_0t}e^{-t^2/2}$
  - Paul小波(mother=1): $\psi(t) = \frac{2^m i^m m!}{\sqrt{\pi(2m)!}}(1-it)^{-(m+1)}$  
  - DOG小波(mother=2): $\psi(t) = \frac{(-1)^{m+1}}{\sqrt{\Gamma(m+1/2)}}D^m[e^{-t^2/2}]$
- 输入参数:
  - Y: 时间序列数据(N点)
  - DT: 时间间隔
  - MOTHER: 母小波类型(0-2)
  - PARAM: 小波参数(k₀、m等)
  - S0,DJ,JTOT: 尺度参数
- 输出参数:
  - WAVE: 小波系数(N×JTOT×2)
  - POWER: 小波功率谱
  - PHASE: 小波相位谱  
  - SIGNIF: 显著性水平
- 功能特点:
  - 基于Torrence & Compo(1998)算法
  - 完整的统计显著性检验
  - 影响锥(COI)计算
  - 全局小波谱(GWS)
  - 适用于时间序列周期分析、气候振荡研究

## 164. wetbulb_profs.f [📝 src](fortran/wetbulb_profs.f)
**功能**: 基于Stipanuk(1973)算法计算湿球温度廓线

**详细说明**:
- 包含多个子程序:
  - WETBULBPROFS: 主接口程序，处理廓线数据
  - TWBPROFS: 核心湿球温度计算
  - OXPROFS,WXPROFS: 位温和混合比计算
  - ESATPROFS: 饱和水汽压计算
  - TSAPROFS,TDAPROFS: 湿绝热线和干绝热线温度
- 输入参数:
  - T: 温度廓线(°C，NPTS点)
  - TD: 露点温度廓线(°C，NPTS点)
  - P: 压力廓线(mb，NPTS点)
- 输出参数:
  - TWB: 湿球温度廓线(°C，NPTS点)
- 计算步骤:
  - 混合比线: $w = \frac{622 \cdot e_s(T_d)}{p - e_s(T_d)}$
  - 干绝热线: $\theta = T \cdot \left(\frac{1000}{p}\right)^{0.286}$
  - 迭代求交点压力: $p_i = p \cdot 2^{0.02(T_{mr}(w,p_i) - T_{da}(\theta,p_i))}$
  - 湿绝热线: $\theta_e = T \cdot \left(\frac{1000}{p}\right)^{0.286} \cdot \exp\left(\frac{2.65 \cdot w}{T}\right)$
- 功能特点:
  - 基于NOAA PROFS系统算法
  - 使用Nordquist(1973)饱和水汽压公式
  - 自动处理缺失值
  - 高精度迭代算法(10次迭代)
  - 适用于大气物理计算、天气预报、探空数据处理

## 165. wgtVertAvg_beta.f [📝 src](fortran/wgtVertAvg_beta.f)
**功能**: 使用Beta因子进行质量加权垂直积分或平均

**详细说明**:
- 包含三个子程序:
  - DWVBETAP1: 一维压力坐标版本
  - DWVBETAP3: 三维压力坐标版本
  - DWVBETAP: 核心Beta因子计算算法
- 输入参数:
  - P: 压力层数组(KLEV或MLON×NLAT×KLEV)
  - X: 待积分变量(MLON×NLAT×KLEV)
  - PSFC: 地面压力(MLON×NLAT)
  - IPUNIT: 压力单位标志(0=mb, 1=Pa)
  - IOPT: 输出选项(0=积分, 1=平均)
  - PTOP,PBOT: 积分上下界
- 输出参数:
  - XVB: 垂直积分或平均结果(MLON×NLAT)
- Beta因子计算:
  - 层厚: $\Delta p(k) = p(k+1) - p(k-1)$
  - Beta因子: $\beta(k) = \frac{\min(p_{bot}, p_{sfc}) - p(k-1)}{p(k+1) - p(k-1)}$
  - 积分: $I = \sum_k X(k) \cdot \beta(k) \cdot \Delta p(k)$
  - 平均: $\bar{X} = \frac{\sum_k X(k) \cdot \beta(k) \cdot \Delta p(k)}{\sum_k \beta(k) \cdot \Delta p(k)}$
- 功能特点:
  - 基于Trenberth(1991)质量守恒诊断方法
  - 自动处理边界层和地形影响
  - 支持部分层积分(PTOP到PBOT)
  - 适用于大气柱总量计算、质量加权平均、垂直积分诊断  
## 166. wk_smooth121.f [📝 src](fortran/wk_smooth121.f)
**功能**: 1-2-1加权平滑滤波器

**详细说明**:
- 包含一个子程序:
  - WKSMOOTH121: 1-2-1权重平滑滤波核心算法
- 输入参数:
  - VV: 输入/输出数据数组(VN)
  - VN: 数组长度
  - NN: 有效数据点数($\leq$VN)
  - SPV: 缺失值标识码
  - DUM: 工作数组(VN)
- 输出参数:
  - VV: 平滑后的数据数组(原地修改)
- 功能特点:
  - 平滑权重公式: $f'(i) = \frac{1 \cdot f(i-1) + 2 \cdot f(i) + 1 \cdot f(i+1)}{4}$
  - 边界处理(第一点): $f'(1) = \frac{3 \cdot f(1) + 1 \cdot f(2)}{4}$
  - 边界处理(最后点): $f'(n) = \frac{1 \cdot f(n-1) + 3 \cdot f(n)}{4}$
  - 保守总和特性: $\sum f'(i) = \sum f(i)$
  - 自动跳过缺失值
  - 适用于Wheeler-Kiladis图分析中的时间序列平滑

## 167. wk_utils.f [📝 src](fortran/wk_utils.f)
**功能**: Wheeler-Kiladis图分析实用工具集

**详细说明**:
- 包含三个子程序:
  - WKSMOOTH121: 1-2-1平滑滤波器(与wk_smooth121.f相同)
  - WKTAPERTOZERO: 时间序列锥化到零值处理
  - WKDETREND: 线性趋势移除
- WKTAPERTOZERO输入参数:
  - TS: 时间序列数据(N)
  - N,NMI,NN,TP: 数组长度、起始索引、有效长度、锥化点数
- WKDETREND输入参数:
  - X2: 时间序列数据(NX,原地修改)
  - NX,VMI,NV: 数组长度、起始位置、有效数据长度
- 功能特点:
  - 锥化窗函数: $w(j) = 0.5 \times (1 - \cos(\frac{(j-1)\pi}{tp}))$, $j = 1,2,\ldots,tp$
  - 锥化公式: $ts'(i) = ts(i) \times w(i)$
  - 线性去趋势: $X'(t) = X(t) - (a + b \cdot t)$
  - 其中 $b = \frac{\sum t \cdot X(t) - n\bar{t}\bar{X}}{\sum t^2 - n\bar{t}^2}$, $a = \bar{X} - b\bar{t}$
  - 满足FFT周期性要求的锥化处理
  - 适用于大气波动频谱分析预处理

## 168. wrf_bint3d.f [📝 src](fortran/wrf_bint3d.f)
**功能**: WRF模式坐标变换和三维双线性插值

**详细说明**:
- 包含多个子程序:
  - DMAPTFORM: WRF地图投影坐标变换
  - DBINT3D: 三维双线性插值主程序
  - DBINT: 二维双线性插值算法
  - DONED: 一维插值核心函数
- DMAPTFORM输入参数:
  - DSKMC,MIYCORS,MJXCORS: 网格间距和维度
  - NPROJ,XLATC,XLONC: 投影类型、中心纬度和经度
  - TRUE1,TRUE2: 真实纬度参数
  - RIY,RJX,RLAT,RLON: 网格坐标和地理坐标
  - IDIR: 转换方向(1=网格→地理，-1=地理→网格)
- DBINT3D输入参数:
  - DATA_IN: 输入三维数据(NX×NY×NZ)
  - OBSII,OBSJJ: 观测点网格坐标(NOBSICRS×NOBSJCRS)
  - ICRS,JCRS: 交叉标志参数
- 输出参数:
  - DATA_OUT: 插值后数据(NOBSICRS×NOBSJCRS×NZ)
- 功能特点:
  - 支持三种地图投影:
    - Lambert Conformal (NPROJ=1): $\text{cone} = \frac{\ln(\cos(\phi_1)/\cos(\phi_2))}{\ln(\tan(\pi/4-\phi_2/2)/\tan(\pi/4-\phi_1/2))}$
    - Polar Stereographic (NPROJ=2): $\text{cone} = 1$
    - Mercator (NPROJ=3): $y = R_e \ln(\tan(\pi/4 + \phi/2))$
  - 四点双线性插值: $f(x,y) = \sum_{i,j} w_{ij} f_{ij}$
  - 权重函数: $w_{ij} = (1-u)(1-v), u(1-v), (1-u)v, uv$
  - 自动处理网格边界和缺失值
  - 适用于WRF模式数据后处理和观测点插值

## 169. wrf_cloud_fracf.f90 [📝 src](fortran/wrf_cloud_fracf.f90)
**功能**: WRF模式云分数计算(Fortran90实现)

**详细说明**:
- 包含两个子程序:
  - DCLOUDFRAC: 固定压力阈值云分数计算(已废弃)
  - DCLOUDFRAC2: 灵活垂直坐标云分数计算(推荐使用)
- DCLOUDFRAC2输入参数:
  - RH: 相对湿度场(EW×NS×NZ，%）
  - VERT: 垂直坐标场(EW×NS×NZ,压力或高度)
  - LOW_THRESH,MID_THRESH,HIGH_THRESH: 低、中、高云阈值
  - VERT_INC_W_HEIGHT: 垂直坐标方向标志(0=递减，1=递增)
  - MSG: 缺失值代码
- 输出参数:
  - LOWC,MIDC,HIGHC: 低云、中云、高云分数(EW×NS)
- 功能特点:
  - 云层分类标准:
    - 低云: $p > 970$ hPa 或 $z < 2$ km
    - 中云: $700 < p < 970$ hPa 或 $2 < z < 6$ km  
    - 高云: $p < 450$ hPa 或 $z > 10$ km
  - 云分数计算:
    - 低云: $CF_{low} = \max(0, \min(1, 4 \times \frac{\text{RH}_{max}}{100} - 3))$
    - 中云: $CF_{mid} = \max(0, \min(1, 4 \times \frac{\text{RH}_{max}}{100} - 3))$
    - 高云: $CF_{high} = \max(0, \min(1, 2.5 \times \frac{\text{RH}_{max}}{100} - 1.5))$
  - 使用OpenMP并行化加速计算
  - 支持任意垂直坐标系统
  - 自动处理地形遮挡情况
  - 适用于WRF模式云量诊断和气候分析

## 170. wrf_constants.f90 [📝 src](fortran/wrf_constants.f90)
**功能**: WRF模式物理常数定义模块

**详细说明**:
- 包含一个Fortran90模块:
  - WRF_CONSTANTS: 完整的物理常数定义集合
- 主要物理常数:
  - 基本常数:
    - $\pi = 3.141592653589793$
    - 地球半径: $R_e = 6.37 \times 10^6$ m
    - 重力加速度: $g = 9.81$ m/s$^2$
  - 热力学常数:
    - 干空气气体常数: $R_d = 287$ J/(kg·K)
    - 水汽气体常数: $R_v = 461.6$ J/(kg·K)
    - 定压比热容: $C_p = 1004.5$ J/(kg·K)
    - 绝热指数: $\gamma = R_d/C_p = 0.286$
  - 水汽参数:
    - 分子量比: $\epsilon = \frac{M_d}{M_v} = 0.622$
    - 饱和水汽压常数: $e_0 = 6.112$ hPa
    - Clausius-Clapeyron常数: $L_1 = 17.67$, $L_2 = 29.65$
  - 云物理参数:
    - 云水吸收系数: $\alpha_w = 0.145$ m$^2$/g
    - 云冰吸收系数: $\alpha_i = 0.272$ m$^2$/g
    - 水密度: $\rho_w = 1000$ kg/m$^3$
- 功能特点:
  - 与WRF模式module_model_constants.F完全兼容
  - 双精度浮点数定义保证数值精度
  - 包含完整的缺失值定义体系
  - 支持多种数据类型的默认填充值
  - 湿绝热过程修正参数:
    - $C_{p,moist} = C_p(1 + 0.887 \times q_v)$
    - $R_{gas,moist} = R_d(1 + 0.608 \times q_v)$
  - 适用于所有WRF相关计算程序的标准化常数引用  
## 171. wgtVertAvg_beta.f [📝 src](fortran/wgtVertAvg_beta.f)
**功能**: 使用Beta因子进行质量加权垂直积分或平均

**详细说明**:
- 包含三个子程序:
  - DWVBETAP1: 一维压力坐标版本
  - DWVBETAP3: 三维压力坐标版本  
  - DWVBETAP: 核心Beta因子计算算法
- 输入参数:
  - P: 压力层数组(KLEV或MLON×NLAT×KLEV)
  - X: 待积分变量(MLON×NLAT×KLEV)
  - PSFC: 地面压力(MLON×NLAT)
  - IPUNIT: 压力单位标志(0=mb, 1=Pa)
  - IOPT: 输出选项(0=积分, 1=平均)
  - PTOP,PBOT: 积分上下界
- 输出参数:
  - XVB: 垂直积分或平均结果(MLON×NLAT)
- 功能特点:
  - Beta因子计算: $\beta(k) = \frac{\min(p_{bot}, p_{sfc}) - p(k-1)}{p(k+1) - p(k-1)}$
  - 层厚: $\Delta p(k) = p(k+1) - p(k-1)$
  - 积分公式: $I = \sum_k X(k) \cdot \beta(k) \cdot \Delta p(k)$
  - 平均公式: $\bar{X} = \frac{\sum_k X(k) \cdot \beta(k) \cdot \Delta p(k)}{\sum_k \beta(k) \cdot \Delta p(k)}$
  - 基于Trenberth(1991)质量守恒诊断方法
  - 自动处理边界层和地形影响
  - 适用于大气柱总量计算、质量加权平均

## 172. wk_smooth121.f [📝 src](fortran/wk_smooth121.f)
**功能**: 1-2-1加权平滑滤波器

**详细说明**:
- 包含一个子程序:
  - WKSMOOTH121: 1-2-1权重平滑滤波核心算法
- 输入参数:
  - VV: 输入/输出数据数组(VN)
  - VN: 数组长度
  - NN: 有效数据点数($\leq$VN)  
  - SPV: 缺失值标识码
  - DUM: 工作数组(VN)
- 输出参数:
  - VV: 平滑后的数据数组(原地修改)
- 功能特点:
  - 平滑权重公式: $f'(i) = \frac{1 \cdot f(i-1) + 2 \cdot f(i) + 1 \cdot f(i+1)}{4}$
  - 边界处理(第一点): $f'(1) = \frac{3 \cdot f(1) + 1 \cdot f(2)}{4}$
  - 边界处理(最后点): $f'(n) = \frac{1 \cdot f(n-1) + 3 \cdot f(n)}{4}$
  - 保守总和特性: $\sum f'(i) = \sum f(i)$
  - 自动跳过缺失值
  - 适用于Wheeler-Kiladis图分析中的时间序列平滑

## 173. wk_utils.f [📝 src](fortran/wk_utils.f)
**功能**: Wheeler-Kiladis图分析实用工具集

**详细说明**:
- 包含三个子程序:
  - WKSMOOTH121: 1-2-1平滑滤波器
  - WKTAPERTOZERO: 时间序列锥化到零值处理
  - WKDETREND: 线性趋势移除
- WKTAPERTOZERO输入参数:
  - TS: 时间序列数据(N)
  - N,NMI,NN,TP: 数组长度、起始索引、有效长度、锥化点数
- WKDETREND输入参数:
  - X2: 时间序列数据(NX,原地修改)
  - NX,VMI,NV: 数组长度、起始位置、有效数据长度
- 功能特点:
  - 锥化窗函数: $w(j) = 0.5 \times (1 - \cos(\frac{(j-1)\pi}{tp}))$, $j = 1,2,\ldots,tp$
  - 锥化公式: $ts'(i) = ts(i) \times w(i)$
  - 线性去趋势: $X'(t) = X(t) - (a + b \cdot t)$
  - 斜率: $b = \frac{\sum t \cdot X(t) - n\bar{t}\bar{X}}{\sum t^2 - n\bar{t}^2}$, $a = \bar{X} - b\bar{t}$
  - 满足FFT周期性要求的锥化处理
  - 适用于大气波动频谱分析预处理

## 174. wrf_bint3d.f [📝 src](fortran/wrf_bint3d.f)
**功能**: WRF模式坐标变换和三维双线性插值

**详细说明**:
- 包含多个子程序:
  - DMAPTFORM: WRF地图投影坐标变换
  - DBINT3D: 三维双线性插值主程序
  - DBINT: 二维双线性插值算法
  - DONED: 一维插值核心函数
- DMAPTFORM输入参数:
  - DSKMC,MIYCORS,MJXCORS: 网格间距和维度
  - NPROJ,XLATC,XLONC: 投影类型、中心纬度和经度
  - TRUE1,TRUE2: 真实纬度参数
  - RIY,RJX,RLAT,RLON: 网格坐标和地理坐标
  - IDIR: 转换方向(1=网格→地理，-1=地理→网格)
- DBINT3D输入参数:
  - DATA_IN: 输入三维数据(NX×NY×NZ)
  - OBSII,OBSJJ: 观测点网格坐标(NOBSICRS×NOBSJCRS)
  - ICRS,JCRS: 交叉标志参数
- 输出参数:
  - DATA_OUT: 插值后数据(NOBSICRS×NOBSJCRS×NZ)
- 功能特点:
  - 支持三种地图投影:
    - Lambert Conformal (NPROJ=1): $\text{cone} = \frac{\ln(\cos(\phi_1)/\cos(\phi_2))}{\ln(\tan(\pi/4-\phi_2/2)/\tan(\pi/4-\phi_1/2))}$
    - Polar Stereographic (NPROJ=2): $\text{cone} = 1$
    - Mercator (NPROJ=3): $y = R_e \ln(\tan(\pi/4 + \phi/2))$
  - 四点双线性插值: $f(x,y) = \sum_{i,j} w_{ij} f_{ij}$
  - 权重函数: $w_{ij} = (1-u)(1-v), u(1-v), (1-u)v, uv$
  - 自动处理网格边界和缺失值
  - 适用于WRF模式数据后处理和观测点插值

## 175. wrf_cloud_fracf.f90 [📝 src](fortran/wrf_cloud_fracf.f90)
**功能**: WRF模式云分数计算(Fortran90实现)

**详细说明**:
- 包含两个子程序:
  - DCLOUDFRAC: 固定压力阈值云分数计算(已废弃)
  - DCLOUDFRAC2: 灵活垂直坐标云分数计算(推荐使用)
- DCLOUDFRAC2输入参数:
  - RH: 相对湿度场(EW×NS×NZ，%)
  - VERT: 垂直坐标场(EW×NS×NZ,压力或高度)
  - LOW_THRESH,MID_THRESH,HIGH_THRESH: 低、中、高云阈值
  - VERT_INC_W_HEIGHT: 垂直坐标方向标志(0=递减，1=递增)
  - MSG: 缺失值代码
- 输出参数:
  - LOWC,MIDC,HIGHC: 低云、中云、高云分数(EW×NS)
- 功能特点:
  - 云层分类标准:
    - 低云: $p > 970$ hPa 或 $z < 2$ km
    - 中云: $700 < p < 970$ hPa 或 $2 < z < 6$ km
    - 高云: $p < 450$ hPa 或 $z > 10$ km
  - 云分数计算:
    - 低云: $CF_{low} = \max(0, \min(1, 4 \times \frac{\text{RH}_{max}}{100} - 3))$
    - 中云: $CF_{mid} = \max(0, \min(1, 4 \times \frac{\text{RH}_{max}}{100} - 3))$
    - 高云: $CF_{high} = \max(0, \min(1, 2.5 \times \frac{\text{RH}_{max}}{100} - 1.5))$
  - 使用OpenMP并行化加速计算
  - 支持任意垂直坐标系统
  - 自动处理地形遮挡情况
  - 适用于WRF模式云量诊断和气候分析

## 176. wrf_constants.f90 [📝 src](fortran/wrf_constants.f90)
**功能**: WRF模式物理常数定义模块

**详细说明**:
- 包含一个Fortran90模块:
  - WRF_CONSTANTS: 完整的物理常数定义集合
- 主要物理常数:
  - 基本常数:
    - $\pi = 3.1415926535897932$
    - 地球半径: $R_e = 6.37 \times 10^6$ m
    - 重力加速度: $g = 9.81$ m/s$^2$
  - 热力学常数:
    - 干空气气体常数: $R_d = 287$ J/(kg·K)
    - 水汽气体常数: $R_v = 461.6$ J/(kg·K)
    - 定压比热容: $C_p = 1004.5$ J/(kg·K)
    - 绝热指数: $\gamma = R_d/C_p = 0.286$
  - 水汽参数:
    - 分子量比: $\epsilon = \frac{M_d}{M_v} = 0.622$
    - 饱和水汽压常数: $e_0 = 6.112$ hPa
    - Clausius-Clapeyron常数: $L_1 = 17.67$, $L_2 = 29.65$
  - 云物理参数:
    - 云水吸收系数: $\alpha_w = 0.145$ m$^2$/g
    - 云冰吸收系数: $\alpha_i = 0.272$ m$^2$/g
    - 水密度: $\rho_w = 1000$ kg/m$^3$
- 功能特点:
  - 与WRF模式module_model_constants.F完全兼容
  - 双精度浮点数定义保证数值精度
  - 包含完整的缺失值定义体系
  - 支持多种数据类型的默认填充值
  - 湿绝热过程修正参数:
    - $C_{p,moist} = C_p(1 + 0.887 \times q_v)$
    - $R_{gas,moist} = R_d(1 + 0.608 \times q_v)$
  - 适用于所有WRF相关计算程序的标准化常数引用

## 177. wrf_fctt.f90 [📝 src](fortran/wrf_fctt.f90)
**功能**: WRF云顶温度计算(基于光学厚度法)

**详细说明**:
- 包含一个子程序:
  - WRFCTTCALC: 云顶温度计算主程序
- 输入参数:
  - PRS: 压力场(ew×ns×nz，Pa)
  - TK: 温度场(ew×ns×nz，K)
  - QCI,QCW: 云冰和云水混合比(ew×ns×nz，g/kg)
  - QVP: 水汽混合比(ew×ns×nz，g/kg)
  - GHT: 位势高度(ew×ns×nz，m)
  - TER: 地形高度(ew×ns，m)
  - HAVEQCI: 是否有云冰数据(0/1)
  - FILL_NOCLOUD: 无云填充选项(0/1)
  - OPT_THRESH: 光学厚度阈值(典型值1.0)
- 输出参数:
  - CTT: 云顶温度(ew×ns，K)
  - PF: 全层压力场(ew×ns×nz，Pa)
- 功能特点:
  - 光学厚度计算:
    - 云水: $\tau_w = \alpha_w \times q_{cw} \times \frac{\Delta p}{g}$
    - 云冰: $\tau_i = \alpha_i \times q_{ci} \times \frac{\Delta p}{g}$
    - 总光学厚度: $\tau = \tau_w + \tau_i$
  - 云顶检测: 从模式顶向下积分至$\tau \geq 1.0$
  - 温度插值: 线性插值获得云顶温度
  - 使用OpenMP并行化加速计算
  - 自动处理地形和边界条件
  - 适用于卫星数据验证和云微物理诊断

## 178. wrf_pvo.f90 [📝 src](fortran/wrf_pvo.f90)
**功能**: WRF位涡和绝对涡度计算

**详细说明**:
- 包含两个子程序:
  - DCOMPUTEABSVORT: 绝对涡度计算
  - DCOMPUTEPV: 位涡计算
- DCOMPUTEABSVORT输入参数:
  - U,V: 风场分量(nxp1×ny×nz, nx×nyp1×nz，m/s)
  - MSFU,MSFV,MSFT: 地图比例因子
  - COR: 科氏参数(nx×ny，s$^{-1}$)
  - DX,DY: 网格间距(m)
- DCOMPUTEPV输入参数:
  - U,V: 风场分量
  - THETA: 位温场(nx×ny×nz，K)
  - PRS: 压力场(nx×ny×nz，Pa)
  - 地图投影和网格参数
- 输出参数:
  - AV: 绝对涡度(nx×ny×nz，×10⁵ s$^{-1}$)
  - PV: 位涡(nx×ny×nz，×10⁻$^2$ PVU)
- 功能特点:
  - 绝对涡度公式: $\zeta_{abs} = \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y} + f$
  - 位涡公式: $PV = -g(\frac{\partial \theta}{\partial p}\zeta_{abs} - \frac{\partial v}{\partial p}\frac{\partial \theta}{\partial x} + \frac{\partial u}{\partial p}\frac{\partial \theta}{\partial y})$
  - 地图投影修正: 使用地图比例因子$m$进行坐标变换
  - 有限差分: 使用中心差分计算偏导数
  - 边界处理: 自动处理网格边界的单边差分
  - 使用OpenMP并行化
  - 适用于动力学诊断和天气分析

## 179. wrf_pw.f90 [📝 src](fortran/wrf_pw.f90)  
**功能**: WRF可降水量计算

**详细说明**:
- 包含一个子程序:
  - DCOMPUTEPW: 大气柱可降水量计算
- 输入参数:
  - P: 压力场(nx×ny×nz，Pa)
  - TV: 虚温场(nx×ny×nz，K)
  - QV: 水汽混合比(nx×ny×nz，kg/kg)
  - HT: 位势高度(nx×ny×nzh，m)，nzh=nz+1
- 输出参数:
  - PW: 可降水量(nx×ny，kg/m$^2$或mm)
- 功能特点:
  - 积分公式: $PW = \sum_{k=1}^{nz} \frac{p(k)}{R_d \times T_v(k)} \times q_v(k) \times \Delta h(k)$
  - 其中$\Delta h(k) = h(k+1) - h(k)$为层厚
  - 密度计算: $\rho = \frac{p}{R_d \times T_v}$
  - 质量积分: $PW = \sum \rho \times q_v \times \Delta h$
  - 单位换算: kg/m$^2$ = mm水柱高度
  - 使用OpenMP并行化
  - 垂直积分从地面到模式顶
  - 适用于大气水汽分析和降水预报

## 180. wrf_relhl.f90 [📝 src](fortran/wrf_relhl.f90)
**功能**: WRF风暴相对螺旋度计算

**详细说明**:
- 包含一个子程序:
  - DCALRELHL: 风暴相对螺旋度计算主程序
- 输入参数:
  - U,V: 风场分量(miy×mjx×mkzh，m/s)
  - GHT: 位势高度(miy×mjx×mkzh，m)
  - TER: 地形高度(miy×mjx，m)
  - LAT: 纬度(miy×mjx，度)
  - TOP: 积分上限高度(m，典型值3000m)
- 输出参数:
  - SREH: 风暴相对螺旋度(miy×mjx，m$^2$/s$^2$)
- 功能特点:
  - 风暴移动速度估算:
    - 平均风速: $\vec{U}_{avg} = \frac{\sum_{k} \vec{U}(k) \Delta h(k)}{\sum_{k} \Delta h(k)}$ (3-10km层平均)
    - 风暴速度: $\vec{C} = 0.75 |\vec{U}_{avg}|$，方向偏右30°(北半球)
  - 螺旋度计算: $SRH = -\sum_{k} [(\vec{U}(k) - \vec{C}) \times (\vec{U}(k) - \vec{U}(k+1))]_z$
  - 物理意义: 测量对流风暴入流环境中的流线涡度
  - 阈值参考:
    - SRH < 100 m$^2$/s$^2$: 截止值
    - SRH = 150-299: 可能产生弱龙卷风的超级单体
    - SRH = 300-499: 有利于强龙卷风发展
    - SRH > 450: 暴力龙卷风
  - 使用OpenMP并行化
  - 自动识别南北半球
  - 适用于强对流天气预报和龙卷风预警  

## 181. wrf_write_wps.f [📝 src](fortran/wrf_write_wps.f)
**功能**: WRF预处理系统中间格式文件写入函数

**详细说明**:
- 包含一个子程序:
  - WRITE_INTERMEDIATE_WPS: 将气象数据写入WPS中间格式文件
- 输入参数:
  - OUTPUT_NAME: 输出文件名前缀
  - FIELDIN: 变量名称(长度$\leq$9字符)
  - UNITSIN: 变量单位(长度$\leq$25字符)  
  - DESCIN: 变量描述(长度$\leq$46字符)
  - HDATEIN: 日期时间字符串(格式YYYY-MM-DD_HH:MM:SS)
  - DATA_SOURCE: 数据源标识(长度$\leq$32字符)
  - XLVL: 垂直层/压力层值
  - IPROJ: 地图投影类型(0-5)
  - START_LOCATION: 起始位置标识(长度$\leq$8字符)
  - STARTLAT,STARTLON: 起始纬度和经度
  - DELTALAT,DELTALON: 纬度和经度间隔
  - XLONC: 投影中心经度  
  - TRUELAT1,TRUELAT2: 投影真实纬度
  - NLATS: 高斯投影纬度数
  - DX,DY: 网格间距(m)
  - NX,NY: 网格维度
  - IS_WIND_EARTH_REL: 风场相对地球坐标标志
  - IFV: WPS格式版本号(默认5)
  - XFCST: 预报时间(小时)
  - EARTH_RADIUS: 地球半径(km)
  - SLAB: 二维数据场(NX×NY)
- 支持的地图投影类型:
  - IPROJ=0: 柱面等距投影(经纬度网格)
  - IPROJ=1: 墨卡托投影
  - IPROJ=3: 兰伯特等角投影  
  - IPROJ=4: 高斯投影
  - IPROJ=5: 极射赤面投影
- 输出文件格式:
  - 文件名格式: `$\text{OUTPUT\_NAME}:\text{YYYY-MM-DD\_HH:MM:SS}$`
  - 二进制大端序(big_endian)格式
  - WPS中间格式记录结构:
    - 记录1: 格式版本号IFV
    - 记录2: 通用信息(日期、变量名、单位、描述等)
    - 记录3: 投影相关参数(根据IPROJ变化)
    - 记录4: 风场坐标系标志
    - 记录5: 数据场SLAB
- 功能特点:
  - 支持文件追加模式,同一时间多个变量可写入同一文件
  - 自动处理不同投影的参数写入
  - 兼容WPS v3.0+格式标准
  - 使用大端序保证跨平台兼容性
  - 适用于WRF数值天气预报模式的预处理数据准备  
## 182. writematrix.f [📝 src](fortran/writematrix.f)
**功能**: 多数据类型矩阵格式化写入文件的I/O函数库

**详细说明**:
- 包含五个子程序:
  - WRITEMATRIXI: 整数矩阵写入函数
  - WRITEMATRIXF: 单精度实数矩阵写入函数
  - WRITEMATRIXD: 双精度实数矩阵写入函数
  - WRITEMATRIXB: 字节型整数矩阵写入函数(INTEGER*1)
  - WRITEMATRIXS: 短整型矩阵写入函数(INTEGER*2)
- 输入参数:
  - FNAME: 输出文件名(如为"*"则输出到标准输出)
  - NROW,NCOL: 矩阵行数和列数
  - X: 输入矩阵数据(NCOL×NROW)
  - FMTX: Fortran格式描述符字符串
  - TITLE: 输出标题字符串
  - TITSP: 标题前空格数
  - IOPT: 输出选项标志
- 输出格式选项:
  - IOPT=0: 仅输出矩阵数据
    - 格式: $\text{FORMAT} = \text{"("} || \text{FMTX} || \text{")"}$
  - IOPT=1: 输出行索引+矩阵数据
    - 格式: $\text{FORMAT} = \text{"(I5,1X,"} || \text{FMTX} || \text{")"}$
- 标题格式化:
  - TITSP>0时: $\text{TITFMT} = \text{"("} || \text{TITSP} || \text{"X,A)"}$
  - TITSP$\leq$0时: $\text{TITFMT} = \text{"(X,A)"}$
- 数据类型支持:
  - INTEGER: 标准整数类型
  - REAL: 单精度浮点数
  - DOUBLE PRECISION: 双精度浮点数  
  - INTEGER*1: 字节型整数(-128到127)
  - INTEGER*2: 短整型整数(-32768到32767)
- 功能特点:
  - 支持文件输出和标准输出两种模式
  - 灵活的格式控制,支持用户自定义Fortran格式符
  - 可选的行索引输出(从0开始编号)
  - 标题输出位置和格式可控制
  - 矩阵按行优先顺序输出: $\text{行i} = [X(1,i), X(2,i), \ldots, X(\text{NCOL},i)]$
  - 自动缓冲区刷新确保输出完整性
  - 支持非标准Fortran扩展数据类型(INTEGER*1/2)
  - 适用于科学计算结果的格式化输出和数据交换  

## 183. wrunave_dp.f [📝 src](fortran/wrunave_dp.f)
**功能**: 加权滑动平均滤波算法(双精度)

**详细说明**:
- 包含两个子程序:
  - DWGTRUNAVE: 加权滑动平均主接口程序
  - DWRUNAVX77: 加权滑动平均核心算法实现
- 输入参数:
  - X: 输入时间序列数据(NPTS，原地修改)
  - NPTS: 数据点数
  - WGT: 权重向量(NWGT)
  - NWGT: 滑动窗口大小(权重点数)
  - KOPT: 边界处理选项
  - XMSG: 缺失值代码
  - WORK: 工作数组(长度$\geq$NPTS+2×NHALF，NHALF=NWGT/2)
  - LWORK: 工作数组长度
- 输出参数:
  - X: 滤波后的时间序列(原地修改)
  - IER: 错误代码(0=成功，<0=错误)
- 边界处理选项:
  - KOPT<0: 循环边界条件
    - 示例(NWGT=3): $X(1) = w_1 X(n) + w_2 X(1) + w_3 X(2)$
    - 示例(NWGT=4): $X(1) = w_1 X(n-1) + w_2 X(n) + w_3 X(1) + w_4 X(2)$
  - KOPT=0: 边界点设为缺失值
    - 不足窗口大小的前后点设为XMSG
    - 示例(NWGT=3): $X(1) = \text{XMSG}$，$X(2) = w_1 X(1) + w_2 X(2) + w_3 X(3)$
  - KOPT>0: 反射(对称)边界条件  
    - 示例(NWGT=3): $X(1) = w_1 X(2) + w_2 X(1) + w_3 X(2)$
    - 示例(NWGT=4): $X(1) = w_1 X(3) + w_2 X(2) + w_3 X(1) + w_4 X(2)$
- 加权平均公式:
  - 权重归一化: $\text{WSUM} = \sum_{i=1}^{\text{NWGT}} w_i$
  - 如果WSUM>1: $\text{WSUM} = \frac{1}{\text{WSUM}}$ (归一化因子)
  - 滤波公式: $X'(n) = \text{WSUM} \times \sum_{k=0}^{\text{NWGT}-1} w_{k+1} \times X(n-\text{NHALF}+k)$
- 缺失值处理:
  - 窗口内任何点为缺失值时，输出设为XMSG
  - 自动跳过缺失值点，不参与滤波计算
- 错误处理:
  - IER=-11: NPTS$\leq$0
  - IER=-12: NWGT>NPTS (窗口大小超过数据长度)
- 功能特点:
  - 双精度浮点运算保证数值精度
  - 支持任意长度的权重向量
  - 三种边界处理策略适应不同应用需求
  - 自动权重归一化确保滤波器增益
  - 内存优化的工作数组设计
  - 原地修改节省内存空间
  - 适用于时间序列平滑、噪声抑制、信号滤波  

## 184. xy1pdf77.f [📝 src](fortran/xy1pdf77.f)
**功能**: 一维概率密度函数估计算法(Fortran77实现)

**详细说明**:
- 包含一个子程序:
  - X1PDF77: 计算一维数据的概率密度函数(直方图)
- 输入参数:
  - NX: 数据点总数
  - X: 输入数据数组(NX)
  - XMSG: 缺失值代码
  - NBX: 直方图区间数
  - NBXP1: 区间边界点数(NBX+1)
  - BINXBND: 区间边界数组(NBXP1)
  - IPCNT: 输出单位标志(0=频次,1=百分比)
- 输出参数:
  - PDF: 概率密度函数值(NBX)
  - IER: 错误代码(0=成功,1=无有效数据)
- 算法流程:
  1. 初始化所有PDF区间为0
  2. 统计有效数据点数: $K_x = \sum_{i=1}^{NX} \mathbf{1}_{X(i) \neq \text{XMSG}}$
  3. 数据装箱处理:
     - 对于区间$j$ (j=1,2,...,NBX): 
     - 计数条件: $\text{BINXBND}(j) \leq X(i) < \text{BINXBND}(j+1)$
     - 特殊处理最后边界: $X(i) = \text{BINXBND}(\text{NBX}+1)$时归入最后区间
  4. 频次计算: $\text{PDF}(j) = \sum_{i=1}^{NX} \mathbf{1}_{X(i) \in \text{区间}j}$
- 输出单位转换:
  - IPCNT=0: 输出原始频次
  - IPCNT=1: 输出百分比 $\text{PDF}(j) = 100 \times \frac{\text{PDF}(j)}{K_x}$
- 性能优化特性:
  - 缺失值检测优化: 当$K_x = NX$时跳过缺失值检查
  - 双层循环结构适应大数据集处理
- 边界处理:
  - 左边界包含: $X(i) \geq \text{BINXBND}(j)$
  - 右边界不包含: $X(i) < \text{BINXBND}(j+1)$ 
  - 最右边界特殊包含: $X(i) = \text{BINXBND}(\text{NBX}+1)$
- 功能特点:
  - 支持任意用户定义的区间边界
  - 自动处理缺失值
  - 可输出频次或百分比两种格式
  - 针对大数据集的性能优化
  - 标准直方图统计算法
  - 适用于数据分布分析、统计建模、数据可视化预处理  

## 185. xy2pdf77.f [📝 src](fortran/xy2pdf77.f)
**功能**: 二维联合概率密度函数估计算法(Fortran77实现)

**详细说明**:
- 包含一个子程序:
  - XY2PDF77: 计算二维数据的联合概率密度函数(二维直方图)
- 输入参数:
  - NXY: 数据对总数
  - X,Y: 输入数据数组(NXY)
  - XMSG,YMSG: X和Y的缺失值代码
  - NBY: Y方向直方图区间数
  - MBX: X方向直方图区间数
  - MBXP1: X方向区间边界点数(MBX+1)
  - NBYP1: Y方向区间边界点数(NBY+1)
  - BINXBND: X方向区间边界数组(MBXP1)
  - BINYBND: Y方向区间边界数组(NBYP1)
  - IPCNT: 输出单位标志(0=频次,1=百分比)
- 输出参数:
  - PDF: 二维概率密度函数矩阵(MBX×NBY)
  - IER: 错误代码(0=成功,1=无有效数据对)
- 算法流程:
  1. 初始化所有PDF区间为0
  2. 统计有效数据对数: $K_{xy} = \sum_{i=1}^{NXY} \mathbf{1}_{X(i) \neq \text{XMSG} \land Y(i) \neq \text{YMSG}}$
  3. 二维数据装箱处理:
     - 对于区间$(m,n)$ (m=1,...,MBX; n=1,...,NBY):
     - X计数条件: $\text{BINXBND}(m) \leq X(i) < \text{BINXBND}(m+1)$
     - Y计数条件: $\text{BINYBND}(n) \leq Y(i) < \text{BINYBND}(n+1)$
  4. 联合频次计算: $\text{PDF}(m,n) = \sum_{i=1}^{NXY} \mathbf{1}_{(X(i),Y(i)) \in \text{区间}(m,n)}$
- 输出单位转换:
  - IPCNT=0: 输出原始频次
  - IPCNT=1: 输出百分比 $\text{PDF}(m,n) = 100 \times \frac{\text{PDF}(m,n)}{K_{xy}}$
- 缺失值处理策略:
  - 严格配对: 只有当$X(i) \neq \text{XMSG}$且$Y(i) \neq \text{YMSG}$时数据对才有效
  - 任一变量缺失则整个数据对被排除
- 性能优化特性:
  - 缺失值检测优化: 当$K_{xy} = NXY$时跳过缺失值检查
  - 三层嵌套循环结构: 外层双循环遍历区间,内层循环遍历数据
- 边界处理:
  - 左边界包含: $X(i) \geq \text{BINXBND}(m)$, $Y(i) \geq \text{BINYBND}(n)$
  - 右边界不包含: $X(i) < \text{BINXBND}(m+1)$, $Y(i) < \text{BINYBND}(n+1)$
- 功能特点:
  - 支持任意用户定义的二维区间网格
  - 严格的缺失值配对处理
  - 可输出频次或百分比两种格式
  - 针对大数据集的性能优化
  - 标准二维直方图统计算法
  - 适用于二维数据分布分析、相关性研究、联合概率建模  

## 186. xy2pdf90.f [📝 src](fortran/xy2pdf90.f)
**功能**: 二维联合概率密度函数估计算法(Fortran90优化实现)

**详细说明**:
- 包含一个子程序:
  - XY2PDF90: 计算二维数据的联合概率密度函数(Fortran90向量化版本)
- 输入参数:
  - NXY: 数据对总数
  - X,Y: 输入数据数组(NXY)
  - XMSG,YMSG: X和Y的缺失值代码
  - NBY: Y方向直方图区间数
  - MBX: X方向直方图区间数
  - MBXP1: X方向区间边界点数(MBX+1)
  - NBYP1: Y方向区间边界点数(NBY+1)
  - BINXBND: X方向区间边界数组(MBXP1)
  - BINYBND: Y方向区间边界数组(NBYP1)
  - IPCNT: 输出单位标志(0=频次,1=百分比)
- 输出参数:
  - PDF: 二维概率密度函数矩阵(MBX×NBY)
  - IER: 错误代码(0=成功,1=无有效数据对)
- 算法流程:
  1. Fortran90数组初始化: $\text{PDF} = 0.0$
  2. 向量化有效数据计数: $K_{xy} = \text{COUNT}(X \neq \text{XMSG} \land Y \neq \text{YMSG})$
  3. 向量化二维装箱处理:
     - 对于区间$(m,n)$ (m=1,...,MBX; n=1,...,NBY):
     - 向量化条件: $\text{LOGICAL MASK} = (X \geq \text{BINXBND}(m)) \land (X < \text{BINXBND}(m+1))$ $\land (Y \geq \text{BINYBND}(n)) \land (Y < \text{BINYBND}(n+1))$
  4. 向量化频次计算: $\text{PDF}(m,n) = \text{COUNT}(\text{LOGICAL MASK})$
- 输出单位转换:
  - IPCNT=0: 输出原始频次
  - IPCNT=1: 向量化百分比转换 $\text{PDF} = 100 \times \frac{\text{PDF}}{K_{xy}}$
- Fortran90优化特性:
  - 内置COUNT函数: 替代显式循环计数,提高执行效率
  - 数组赋值语法: $\text{PDF} = 0.0$和$\text{PDF} = \text{UNIT} \times (\text{PDF}/K_{xy})$
  - 逻辑数组操作: 使用.AND.和.NE.等逻辑运算符的向量化
  - 避免显式内层循环: 使用COUNT内置函数替代三层嵌套循环
- 性能改进:
  - 相比Fortran77版本减少循环层数
  - 利用编译器向量化优化
  - 内存访问模式优化
- 缺失值处理策略:
  - 严格配对: 只有当$X(i) \neq \text{XMSG}$且$Y(i) \neq \text{YMSG}$时数据对才有效
  - 向量化缺失值检测和过滤
- 边界处理:
  - 左边界包含: $X(i) \geq \text{BINXBND}(m)$, $Y(i) \geq \text{BINYBND}(n)$
  - 右边界不包含: $X(i) < \text{BINXBND}(m+1)$, $Y(i) < \text{BINYBND}(n+1)$
- 功能特点:
  - Fortran90现代语法和内置函数
  - 向量化操作显著提高性能
  - 简洁的代码结构和更高的可读性
  - 保持与Fortran77版本完全一致的算法逻辑
  - 自动编译器优化支持
  - 适用于大规模二维数据分布分析、高性能统计计算  

## 187. z2geouv_dp.f [📝 src](fortran/z2geouv_dp.f)
**功能**: 位势高度场到地转风场转换算法(双精度)

**详细说明**:
- 包含三个子程序:
  - Z2GEOUV: 主接口程序，自动判断纬度排列顺序
  - ZUVNEW: 北到南数据重排序处理程序  
  - Z2GUV: 地转风核心计算程序(要求南到北排列)
- 输入参数:
  - Z: 位势高度场(MLON×NLAT，单位：gpm)
  - MLON,NLAT: 经纬度网格点数
  - ZMSG: 缺失值代码
  - GLON: 经度数组(MLON)
  - GLAT: 纬度数组(NLAT)
  - IOPT: 边界条件选项(0=非周期性，1=周期性)
- 输出参数:
  - UG,VG: 地转风U和V分量(MLON×NLAT，单位：m/s)
- 地转风基本方程:
  - 地转平衡: $f \vec{V_g} = -\nabla \Phi \times \hat{k}$
  - U分量: $U_g = -\frac{g}{f} \frac{\partial Z}{\partial y} = -\frac{g}{f} \frac{1}{R} \frac{\partial Z}{\partial \phi}$
  - V分量: $V_g = \frac{g}{f} \frac{\partial Z}{\partial x} = \frac{g}{f} \frac{1}{R \cos\phi} \frac{\partial Z}{\partial \lambda}$
- 物理常数和参数:
  - 重力加速度: $g = 9.80616 \text{ m/s}^2$
  - 地球半径: $R_e = 6371220 \text{ m}$
  - 地球自转角速度: $\Omega = 7.292 \times 10^{-5} \text{ rad/s}$
  - 科里奥利参数: $f = 2\Omega \sin\phi$
- 计算方法:
  - 重力/科氏力比值: $\frac{g}{f} = \frac{g}{2\Omega \sin\phi}$
  - 空间差分距离: 
    - 经向: $\Delta x = 2R \Delta\lambda \cos\phi$ (单位：m)
    - 纬向: $\Delta y = R \Delta\phi$ (单位：m)
  - 中心差分格式: $\frac{\partial Z}{\partial x} \approx \frac{Z(i+1,j) - Z(i-1,j)}{2\Delta x}$
- 边界处理策略:
  - 纬度边界: 使用单侧差分，距离权重减半
  - 经度边界: 
    - IOPT=1: 周期性边界条件 $Z(0,j) = Z(\text{MLON},j)$
    - IOPT=0: 复制相邻点值，V分量乘以2补偿单侧差分
- 特殊区域处理:
  - 赤道附近: 当$|\phi| < 10^{-5}$时，科氏力项设为缺失值
  - 极地附近: 当$|\phi| > 89.999°$时，经向距离设为0
- 数据排列要求:
  - 算法要求纬度从南到北排列
  - 自动检测输入数据排列顺序
  - 北到南数据自动重排序处理
- 功能特点:
  - 双精度浮点计算保证精度
  - 自动数据排列顺序检测和处理
  - 灵活的边界条件设置
  - 完整的地球物理参数集
  - 准确的球面几何计算
  - 适用于数值天气预报、大气动力学分析、风场诊断

## 188. zon_mpsi.f [📝 src](fortran/zon_mpsi.f)
**功能**: 经向平均流函数计算算法(大气环流分析)

**详细说明**:
- 包含两个子程序:
  - DZPSIDRV: 主驱动程序，分配工作数组
  - DZONMPSI: 经向平均流函数核心计算程序
- 输入参数:
  - V: 经向风场(MLON×NLAT×KLEV，单位：m/s)
  - LAT: 纬度数组(NLAT，单位：度)
  - P: 压力层数组(KLEV，单位：Pa)
  - PS: 地面气压(MLON×NLAT，单位：Pa)  
  - VMSG: 缺失值代码
  - MLON,NLAT,KLEV: 经度、纬度、垂直层数
- 输出参数:
  - ZMPSI: 经向平均流函数(NLAT×KLEV，单位：kg·m/s$^2$)
- 流函数计算基本方程:
  - 连续性方程: $\frac{\partial \bar{v}}{\partial p} + \frac{\partial \bar{\omega}}{\partial y} = 0$
  - 流函数定义: $\bar{v} = -\frac{\partial \Psi}{\partial p}$，$\bar{\omega} = \frac{1}{a\cos\phi}\frac{\partial \Psi}{\partial \phi}$
  - 积分关系: $\Psi(\phi,p) = \int_{p_{top}}^{p} \bar{v}(\phi,p') dp'$
- 计算方法和参数:
  - 地球物理常数:
    - 重力加速度: $g = 9.80616 \text{ m/s}^2$
    - 地球半径: $a = 6.37122 \times 10^6 \text{ m}$
    - 转换常数: $C = \frac{2\pi a}{g}$
  - 纬度权重: $C(\phi) = C \times \cos\phi$
- 算法流程:
  1. 压力层扩展: 创建半层压力 $P_{half}(k) = \frac{P(k)+P(k+1)}{2}$
  2. 压力差计算: $\Delta P(k) = P(k+1) - P(k-1)$
  3. 经向风场预处理: 将$P > P_s$的格点设为缺失值
  4. 纬向平均计算: $\bar{v}(\phi,p) = \frac{1}{N_{lon}} \sum_{i=1}^{N_{lon}} v(i,\phi,p)$
  5. 垂直积分: $\Psi(\phi,p) = -C(\phi) \sum_{k=1}^{p} \bar{v}(\phi,k) \Delta P(k)$
- 边界条件处理:
  - 顶层边界: $P_{top} = 500 \text{ Pa}$ (平流层下边界)
  - 底层边界: $P_{bot} = 100500 \text{ Pa}$ (地面以下)
  - 下边界条件: $\Psi(\phi,P_{bot}) = 0$ (质量守恒要求)
- 流函数物理意义:
  - 正值: 顺时针环流(北半球)
  - 负值: 逆时针环流(北半球)
  - 等值线: 代表平均经向环流的流线
- 数值处理特性:
  - 自动缺失值传播: 任一层级缺失值会影响积分结果
  - 质量守恒校正: 强制底边界流函数为零
  - CSM约定: 输出值取负号以符合气候模式约定
- 功能特点:
  - 基于质量连续性方程的严格推导
  - 支持任意垂直坐标系统
  - 自动处理地形以下区域
  - 完整的大气环流诊断
  - 适用于Hadley环流、Walker环流分析
  - 支持气候模式输出的标准化处理
  - 适用于大气环流诊断、气候变化研究、经向热量输送分析

## 189. zregr.f [📝 src](fortran/zregr.f)
**功能**: 多元线性回归分析算法(含标准化系数)

**详细说明**:
- 包含两个子程序:
  - DZREGR1: 多元回归主程序，计算标准化系数
  - DZREGR2: 回归核心算法，最小二乘法求解
- 输入参数:
  - Y: 因变量时间序列(N)
  - X: 自变量矩阵(N×M)
  - N: 观测样本数
  - M: 自变量个数
  - XMSG,YMSG: X和Y的缺失值代码
  - 工作数组: WK,YY,COV,XSD,XMEAN,A,AINV,S
- 输出参数:
  - C: 原始回归系数(M)
  - CNORM: 标准化回归系数(M) 
  - CON: 回归常数项
  - RESID: 残差序列(N)
- 回归模型方程:
  - 原始模型: $Y(t) = C_0 + \sum_{j=1}^{M} C_j X_j(t) + \varepsilon(t)$
  - 标准化模型: $\frac{Y-\bar{Y}}{\sigma_Y} = \sum_{j=1}^{M} C_{norm,j} \frac{X_j-\bar{X_j}}{\sigma_{X_j}}$
- 最小二乘法求解:
  - 矩阵形式: $\vec{Y} = \mathbf{X}\vec{C} + \vec{\varepsilon}$
  - 正规方程: $\mathbf{A}\vec{C} = \vec{S}$ 其中 $\mathbf{A} = \mathbf{X}^T\mathbf{X}$，$\vec{S} = \mathbf{X}^T\vec{Y}$
  - 解向量: $\vec{C} = \mathbf{A}^{-1}\vec{S} = (\mathbf{X}^T\mathbf{X})^{-1}\mathbf{X}^T\vec{Y}$
- 统计量计算:
  - 常数项: $C_0 = \bar{Y} - \sum_{j=1}^{M} C_j \bar{X_j}$
  - 标准化系数: $C_{norm,j} = C_j \frac{\sigma_{X_j}}{\sigma_Y}$
  - 残差: $\text{RESID}(i) = Y(i) - \hat{Y}(i)$
  - 残差方差: $\sigma^2 = \frac{1}{N-M} \sum_{i=1}^{N} \text{RESID}(i)^2$
- 协方差矩阵计算:
  - 系数协方差: $\text{COV}(\vec{C}) = \sigma^2 (\mathbf{X}^T\mathbf{X})^{-1}$
  - 标准误差: $\text{SE}(C_j) = \sqrt{\text{COV}(j,j)}$
- 数值方法特性:
  - 高斯消元法求解矩阵逆
  - 主元选择避免数值奇异性
  - 自动缺失值处理和传播
- 缺失值处理策略:
  - 严格逐行处理: 任一变量缺失则整行数据排除
  - 有效样本计数: 仅使用完整数据进行回归
  - 缺失值标记传播到残差序列
- 统计诊断功能:
  - 决定系数: $R^2 = 1 - \frac{\sum \text{RESID}^2}{\sum (Y-\bar{Y})^2}$
  - F统计量和显著性检验支持
  - 残差分析和模型验证
- 标准化回归意义:
  - 消除量纲影响，便于比较各变量相对重要性
  - 标准化系数绝对值大小反映变量贡献度
  - 适用于多变量重要性排序分析
- 功能特点:
  - 完整的多元线性回归实现
  - 同时提供原始和标准化系数
  - 严格的数值稳定性保证
  - 全面的统计诊断输出
  - 灵活的缺失值处理机制
  - 适用于气候指数建模、预报方程建立、多变量相关分析

## 190. ompgen.F90
**功能**: OpenMP并行代码自动生成器(Fortran90)

**详细说明**:
- 功能特点:
  - OpenMP指令自动插入
  - Fortran90自由格式源码处理
  - 并行代码模式生成
  - 支持模板化代码生成
  - 优化循环并行化

## 191. triple2grid.f [📝 src](fortran/triple2grid.f)
**功能**: 三元组散点数据到规则网格插值算法

**详细说明**:
- 包含三个子程序:
  - TRIPLE2GRID1: 主插值接口程序
  - TRIP2GRD2: 快速最近邻插值算法
  - TRIP2GRD3: 全搜索最近邻插值算法
- 输入参数:
  - KZ: 输入数据点总数
  - XI,YI,ZI: 输入三元组数据(x坐标,y坐标,z值)
  - ZMSG: 缺失值代码
  - MX,NY: 输出网格维度
  - GX,GY: 输出网格坐标数组
  - DOMAIN: 域扩展系数
  - LOOP: 算法选择标志(0=TRIP2GRD2, 1=TRIP2GRD3)
  - METHOD: 距离计算方法(0=欧几里得距离, 1=球面距离)
  - DISTMX: 最大搜索距离
- 输出参数:
  - GRID: 插值后的规则网格(MX×NY)
  - IER: 错误代码
- 插值算法原理:
  - 欧几里得距离: $d = \sqrt{(x_i - x_j)^2 + (y_i - y_j)^2}$
  - 球面距离: $d = R_e \times \arccos(\sin\phi_i \sin\phi_j + \cos\phi_i \cos\phi_j \cos(\lambda_i - \lambda_j))$
  - 最近邻分配: 将每个观测点分配给距离最近的网格点
- 算法优化特性:
  - 等间距网格检测: 自动检测X/Y方向是否等间距以优化搜索
  - 缺失值预处理: 自动剔除缺失值点
  - 外点处理: 通过域扩展处理网格边界外的数据点
  - 精确匹配优先: 坐标完全匹配的点优先处理
- 距离计算选项:
  - METHOD=0: 笛卡尔坐标系，适用于小区域投影数据
  - METHOD=1: 球面坐标系，适用于经纬度数据，使用地球半径$R_e = 6371.22$ km
- 边界处理策略:
  - 域内数据: 直接在目标网格上插值
  - 域外数据: 创建扩展网格，然后提取目标区域
  - 扩展因子: DOMAIN参数控制边界扩展范围
- 性能对比:
  - TRIP2GRD2: 快速算法，适用于大多数情况
  - TRIP2GRD3: 全搜索算法，计算量大但覆盖性好
- 功能特点:
  - 支持不规则分布的散点数据
  - 自适应网格间距检测
  - 双精度浮点计算精度
  - 灵活的距离计算模式
  - 智能的边界外点处理
  - 适用于气象站观测数据网格化、海洋数据插值、地理信息系统

## 192. grid2triple.f [📝 src](fortran/grid2triple.f)
**功能**: 规则网格数据转换为三元组(三列)格式

**详细说明**:
- 包含一个子程序:
  - GRID2TRIPLE: 将2D网格数据转换为三元组列表
- 输入参数:
  - X(MX): X方向坐标数组
  - Y(NY): Y方向坐标数组  
  - Z(MX,NY): 2D数据网格
  - MX,NY: 网格维度
  - LDMAX: 输出数组最大容量(通常为MX×NY)
  - ZMSG: 缺失值代码
- 输出参数:
  - D(LDMAX,3): 三元组数据数组，每行包含[x,y,z]值
  - LD: 实际有效数据点个数(排除缺失值)
  - IER: 错误代码(-10表示所有值都是缺失值)
- 转换过程:
  - 按行优先顺序遍历2D网格: Z(1,1), Z(2,1), ..., Z(MX,1), Z(1,2), ...
  - 跳过缺失值点: 仅保存 Z(m,n) $\neq$ ZMSG 的数据点
  - 输出格式: D(i,1)=X(m), D(i,2)=Y(n), D(i,3)=Z(m,n)
- 数据结构转换:
  - 输入: 结构化网格 Z[MX×NY] + 坐标向量 X[MX], Y[NY]
  - 输出: 非结构化点云 [(x₁,y₁,z₁), (x₂,y₂,z₂), ..., (xₙ,yₙ,zₙ)]
- 应用场景:
  - 网格数据预处理: 为散点插值算法准备输入数据
  - 数据格式转换: 从规则网格转为点云格式
  - 统计分析: 便于对网格数据进行点对点统计
  - 可视化: 转换为适合散点图的数据格式
- 功能特点:
  - 自动过滤缺失值，提高数据质量
  - 保持坐标-数值对应关系
  - 内存高效的单遍扫描算法
  - 支持任意大小的规则网格
  - 适用于气候数据分析、地理信息处理、科学可视化

## 193. linint2.f [📝 src](fortran/linint2.f)
**功能**: 一维和二维线性插值算法集合

**详细说明**:
- 包含多个子程序:
  - DLININT1: 一维分段线性插值主程序
  - DLININT2: 二维双线性插值主程序  
  - DLININT2PTS: 二维点对点插值程序
  - DLIN2INT1: 一维插值核心算法
  - DLINT2XY: 二维插值核心算法
  - DLINCYC: 循环边界处理
  - DMONOID1/DMONOID2: 单调性检验
  - ESTFOW: 缺失值估算算法
- 一维插值 (DLININT1):
  - 输入参数: XI(NXI), FI(NXI), XO(NXO), ICYCX, IOPT
  - 输出参数: FO(NXO), IER
  - 支持循环边界: ICYCX=0(非循环), $\neq$0(循环)
  - 插值选项: IOPT=0(保持缺失区域), IOPT=1(尽可能填充)
- 二维插值 (DLININT2):
  - 输入参数: XI(NXI), YI(NYI), FI(NXI,NYI), XO(NXO), YO(NYO)
  - 输出参数: FO(NXO,NYO), IER
  - 算法流程: 先X方向插值 → 再Y方向插值
  - 双线性插值公式: $f(x,y) = f(x_1,y_1)(1-t_x)(1-t_y) + f(x_2,y_1)t_x(1-t_y) + f(x_1,y_2)(1-t_x)t_y + f(x_2,y_2)t_x t_y$
  - 其中: $t_x = \frac{x-x_1}{x_2-x_1}$, $t_y = \frac{y-y_1}{y_2-y_1}$
- 点对点插值 (DLININT2PTS):
  - 输入坐标对: XO(NXYO), YO(NXYO)
  - 输出: FO(NXYO) 
  - 适用于不规则分布的目标点
- 循环边界处理特性:
  - 自动扩展边界: 在数组两端添加虚拟点
  - 保持周期性: 左边界=右边界值，右边界=左边界值
  - 间距估算: 根据现有网格间距推算边界间距
- 缺失值处理策略:
  - 精确匹配: 坐标完全相同时直接赋值
  - 完整插值: 四个角点都有效时执行标准双线性插值
  - 部分插值: 部分角点缺失时使用距离加权平均
  - 距离权重: $w_{i,j} = \frac{1}{\sqrt{(x_i-x_0)^2+(y_j-y_0)^2}}$
- 单调性验证:
  - DMONOID1: 检验单序列单调性(递增/递减)
  - DMONOID2: 检验两序列单调性一致性
  - 支持严格单调递增和严格单调递减
- 功能特点:
  - 高效的分段线性插值算法
  - 完整的循环边界支持
  - 智能的缺失值处理
  - 严格的输入验证和错误检查
  - 双精度浮点计算精度
  - 适用于气象数据插值、海洋学网格化、地球科学数值模拟

## 194. linmsg_dp.f [📝 src](fortran/linmsg_dp.f)
**功能**: 时间序列缺失值线性插值填充算法

**详细说明**:
- 包含一个子程序:
  - DLINMSG: 一维时间序列缺失值线性插值主程序
- 输入参数:
  - X(NPTS): 输入时间序列数组(可能包含缺失值)
  - NPTS: 序列长度
  - XMSG: 缺失值代码
  - MFLAG: 端点处理标志
    - MFLAG < 0: 开始和结束的缺失值设为最近非缺失值
    - MFLAG $\geq$ 0: 开始和结束的缺失值保持为缺失值
  - MPTCRT: 最大连续缺失值阈值
    - 当连续缺失值超过此阈值时不进行插值
    - 通常设为NPTS以尽可能多地插值
- 输出结果:
  - X(NPTS): 插值后的时间序列(原地修改)
  - 返回状态码通过内部变量跟踪
- 算法流程:
  1. **全缺失检查**: 如果整个序列都是缺失值，直接返回
  2. **缺失段识别**: 扫描序列标记连续缺失值段的起始(NSTRT)和结束(NEND)位置
  3. **阈值检查**: 如果缺失段长度 > MPTCRT，跳过该段不插值
  4. **开始端插值**: 如果序列开始就是缺失值(NSTRT=1)
     - MFLAG < 0: 用第一个有效值填充
     - MFLAG $\geq$ 0: 保持缺失状态
  5. **中间段插值**: 使用线性插值公式
     - 斜率: $\text{slope} = \frac{X(N) - X(\text{NBASE})}{N - \text{NBASE}}$
     - 插值: $X(NN) = X(\text{NBASE}) + \text{slope} \times (NN - \text{NBASE})$
  6. **结束端插值**: 如果序列结束是缺失值(NEND=NPTS)
     - MFLAG < 0: 用最后一个有效值填充
     - MFLAG $\geq$ 0: 保持缺失状态
- 插值特性:
  - 保持数据的时间连续性
  - 线性趋势保持: 插值点严格位于前后有效点连线上
  - 边界处理灵活: 可选择填充或保持缺失状态
  - 阈值控制: 避免跨越过长缺失段的不合理插值
- 应用场景:
  - 气象观测数据预处理: 填补仪器故障期间的数据空隙
  - 时间序列分析: 为统计分析准备连续数据
  - 数据质量控制: 合理范围内的缺失值修复
  - 模型输入准备: 为数值模型提供完整时间序列
- 功能特点:
  - 原地操作，内存使用高效
  - 支持任意长度的时间序列
  - 灵活的端点处理策略
  - 可配置的插值阈值控制
  - 双精度浮点计算精度
  - 适用于各种时序数据的预处理和质量控制

## 195. moc_loops.f [📝 src](fortran/moc_loops.f)
**功能**: 海洋经向翻转环流(MOC)区域积分计算

**详细说明**:
- 包含一个子程序:
  - MOCLOOPS: 海洋经向翻转环流的纬度带和深度层积分主程序
- 输入参数:
  - NYAUX: 辅助纬度网格点数
  - MLON,NLAT,KDEP: 经度、纬度、深度维度
  - NRX: 区域数量(1=全球，2=大西洋)
  - TLAT(MLON,NLAT): 每个网格点的实际纬度值
  - LAT_AUX_GRID(NYAUX): 辅助纬度网格边界
  - RMLAK(MLON,NLAT,2): 区域掩膜数组
    - RMLAK(:,:,1)=1: 全球海洋点
    - RMLAK(:,:,2)=1: 大西洋海洋点
  - WORK1,WORK2,WORK3(MLON,NLAT,KDEP): 三个输入的三维海洋数据场
  - WMSG: 缺失值代码
- 输出参数:
  - TMP1,TMP2,TMP3(NYAUX,KDEP,2): 纬度-深度-区域积分结果
    - 第一维: 纬度带索引
    - 第二维: 深度层索引  
    - 第三维: 区域索引(1=全球，2=大西洋)
- 算法流程:
  1. **初始化**: 将所有输出数组TMP1,TMP2,TMP3设置为0
  2. **全球积分** (NR=1):
     - 遍历纬度带: LAT_AUX_GRID(NY-1) $\leq$ TLAT < LAT_AUX_GRID(NY)
     - 仅处理全球海洋点: RMLAK(ML,NL,1)=1
     - 跳过缺失值: WORK1(ML,NL,KD) $\neq$ WMSG
     - 累积求和: TMP1(NY,KD,1) += WORK1(ML,NL,KD)
  3. **大西洋积分** (NR=2):
     - 相同的纬度带分组逻辑
     - 仅处理大西洋海洋点: RMLAK(ML,NL,2)=1
     - 累积结果到TMP1(NY,KD,2), TMP2(NY,KD,2), TMP3(NY,KD,2)
- MOC计算原理:
  - 经向翻转环流反映海洋的南北向质量输送
  - 按纬度带积分消除了经向(东西向)变化
  - 按深度层分离了不同水深的环流结构
  - 区域分离可对比全球和局部海盆的环流差异
- 数据处理特性:
  - 纬度带自动分组: 根据实际纬度值和辅助网格自动归类
  - 区域选择性积分: 通过掩膜数组实现海陆分离和海盆分离
  - 缺失值安全处理: 自动跳过无效数据点
  - 三场并行计算: 可同时处理速度、温度、盐度等多个物理量
- 应用场景:
  - 气候模式后处理: 计算海洋模式的经向翻转环流强度
  - 海洋环流分析: 研究不同纬度和深度的海水交换
  - 气候变化监测: 追踪长期环流变化趋势  
  - 模式对比验证: 评估不同海洋模式的环流表现
- 功能特点:
  - 高效的三重嵌套循环算法
  - 灵活的区域掩膜支持
  - 同时处理多个物理场
  - 自动纬度带分组机制
  - 双精度浮点计算精度
  - 适用于全球海洋模式的MOC诊断分析

## 196. rcm2points.f [📝 src](fortran/rcm2points.f)
**功能**: 区域气候模式到散点的插值算法

**详细说明**:
- 包含一个子程序:
  - DRCM2POINTS: 将规则/非规则2D网格数据插值到指定散点位置
- 输入参数:
  - NGRD: 数据场层数(可同时处理多个变量)
  - NXI,NYI: 输入网格维度
  - XI(NXI,NYI),YI(NXI,NYI): 输入网格的经纬度坐标(可为非规则网格)
  - FI(NXI,NYI,NGRD): 多层输入数据场
  - NXYO: 输出散点个数
  - XO(NXYO),YO(NXYO): 目标散点的经纬度坐标
  - XMSG: 缺失值代码
  - OPT: 插值方法选项(0/1=距离反比插值, 2=双线性插值)
  - NCRIT: 有效点阈值(通常$\geq$3)
  - KVAL: 网格步长参数(通常=1)
- 输出参数:
  - FO(NXYO,NGRD): 插值后的多层散点数据
  - IER: 错误代码
- 三阶段插值算法:
  1. **精确匹配阶段**: 
     - 寻找坐标完全相同的点: XO = XI 且 YO = YI
     - 直接赋值: FO(NXY,NG) = FI(IX,IY,NG)
  2. **局部插值阶段**:
     - 在2×2网格单元内寻找目标点
     - 条件: XI(IX,IY) $\leq$ XO $\leq$ XI(IX+K,IY) 且 YI(IX,IY) $\leq$ YO $\leq$ YI(IX,IY+K)
     - 双线性插值(OPT=2): $w_{i,j} = (1-t_x)(1-t_y), t_x(1-t_y), (1-t_x)t_y, t_x t_y$
     - 距离反比插值(OPT$\neq$2): $w_{i,j} = \frac{1}{d_{i,j}^2}$，其中$d_{i,j}$为球面距离
  3. **全局搜索阶段**:
     - 对仍未插值的点执行半径搜索
     - 搜索半径: DKM = 5°纬度 $\approx$ 555公里
     - 条件: |YI - YO| $\leq$ 5° 且球面距离 $\leq$ 555km
     - 距离反比加权: $w = \frac{1}{d^2}$，$FO = \frac{\sum w_i \cdot FI_i}{\sum w_i}$
- 距离计算:
  - 使用DGCDIST函数计算球面大圆距离
  - 公式: $d = R_e \times \arccos(\sin\phi_1 \sin\phi_2 + \cos\phi_1 \cos\phi_2 \cos(\lambda_1-\lambda_2))$
  - 地球半径: RE = 6371 km
- 质量控制特性:
  - 单调性检验: 验证输入网格的经纬度单调性
  - 有效点检验: 要求至少NCRIT个有效邻点才进行插值
  - 缺失值处理: 自动跳过缺失值并调整权重
  - 多层同步: 所有数据层使用相同的插值权重
- 应用场景:
  - 区域气候模式后处理: 将模式网格数据插值到观测站点
  - 模式-观测对比: 从模式网格提取观测点位置的数值
  - 数据同化: 为数据同化系统准备观测算子
  - 敏感性分析: 提取特定位置的模式输出
- 功能特点:
  - 支持非规则网格输入(曲线网格)
  - 多种插值方法可选(双线性/距离反比)
  - 三阶段渐进式插值策略
  - 同时处理多个数据场
  - 严格的质量控制和错误检查
  - 适用于区域气候模式、海洋模式、大气模式的散点插值

## 197. rcm2rgrid.f [📝 src](fortran/rcm2rgrid.f)
**功能**: 区域气候模式与规则网格双向插值算法

**详细说明**:
- 包含三个子程序:
  - DRCM2RGRID: 曲线网格到规则网格的插值
  - DRGRID2RCM: 规则网格到曲线网格的插值  
  - DGCDIST: 球面大圆距离计算函数
- DRCM2RGRID 插值特性:
  - 输入: 曲线网格 XI(NXI,NYI), YI(NXI,NYI), FI(NXI,NYI,NGRD)
  - 输出: 规则网格 XO(NXO), YO(NYO), FO(NXO,NYO,NGRD)
  - 三阶段插值策略:
    1. **近似精确匹配**: 容差EPS=1×10⁻$^4$，寻找几乎重叠的点
    2. **距离反比插值**: 在2×2曲线网格单元内使用球面距离反比权重
    3. **线性插值填充**: 对剩余缺失点沿经度方向执行线性插值
- DRGRID2RCM 插值特性:
  - 输入: 规则网格 XI(NXI), YI(NYI), FI(NXI,NYI,NGRD)
  - 输出: 曲线网格 XO(NXO,NYO), YO(NXO,NYO), FO(NXO,NYO,NGRD)
  - 双阶段插值策略:
    1. **近似精确匹配**: 容差EPS=1×10⁻$^3$
    2. **双线性/距离反比插值**: 
       - 四个角点完整: 标准双线性插值 FBLI()
       - 部分缺失: 距离反比权重插值
- 双线性插值公式:
  - 内联函数: FLI(Z1,Z2,SLOPE) = Z1 + SLOPE×(Z2-Z1)
  - 二维扩展: FBLI(Z1,Z2,Z3,Z4,SLPX,SLPY) = 先X插值，再Y插值
  - 权重计算: SLPX = (XO-XI₁)/(XI₂-XI₁), SLPY = (YO-YI₁)/(YI₂-YI₁)
- DGCDIST 球面距离算法:
  - 输入: 两点经纬度 (RLAT1,RLON1), (RLAT2,RLON2)
  - 输出单位选择: IU=1(弧度), IU=2(度), IU=3(米), IU=4(公里)
  - 计算公式: Haversine变种，处理了经度跨180°的情况
  - 特殊优化: 同点检测直接返回0距离
  - 地球半径: 6371.22 km
- 插值算法优化:
  - 单调性预检: 确保输入输出网格的坐标单调性
  - 精确匹配优先: 避免不必要的数值插值误差
  - 阈值控制: NCRIT参数控制最少有效邻点数
  - 缺失值安全: 完整的缺失值检查和权重归一化
- 曲线网格处理特性:
  - 支持任意曲线坐标系(如Lambert投影、球面网格)
  - 自动检测网格单元包含关系
  - 处理网格畸变和不规则间距
  - 三阶段渐进填充策略确保高覆盖率
- 应用场景:
  - 区域气候模式后处理: RCM输出网格化到标准纬度-经度网格
  - 模式嵌套: 全球模式数据插值到区域模式网格
  - 观测同化: 卫星观测数据重新网格化
  - 数据标准化: 不同模式间的网格统一
- 功能特点:
  - 双向插值能力(曲线↔规则)
  - 多种插值方法自适应选择
  - 高精度球面几何计算
  - 批量多层数据处理
  - 鲁棒的质量控制机制
  - 适用于数值天气预报、气候建模、海洋学研究


# 总结

本文档涵盖了NCL项目中的197个Fortran源文件，包括：
- 173个Fortran 77格式文件(.f)
- 21个Fortran 90/95格式文件(.f90)
- 1个Fortran 90模块文件(.F90)
- 1个Python接口定义文件(.pyf)

这些文件涉及气象学、海洋学、气候学、数值分析、统计学等多个领域的科学计算算法。

这些文件涵盖了以下主要功能领域：

## 核心功能分类：

### 1. 气象学与大气科学 (30%)
- 温度、湿度、压力转换
- 大气物理过程计算
- WRF模式专用函数
- 垂直坐标转换
- 大气动力学计算

### 2. 数值分析与计算 (25%)
- FFT变换和频域分析
- 插值和外推算法
- 数值积分与微分
- 矩阵运算和线性代数
- 求解器和优化算法

### 3. 统计分析 (20%)
- 概率分布函数
- 回归分析和相关性
- 时间序列分析
- 假设检验
- 聚类和分类

### 4. 数据处理与格点化 (15%)
- 网格转换和重映射
- 缺失值处理
- 数据平均和聚合
- 质量控制
- 文件I/O操作

### 5. 球面几何与地球科学 (10%)
- 球谐分析
- 地理坐标转换
- 大圆计算
- 地球物理参数

