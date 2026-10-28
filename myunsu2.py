import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import itertools
import jenkspy  # pip install jenkspy
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

# ==========================================
# 0. 环境设置与路径
# ==========================================
output_dir = 'Results_Great/Hazard_Unsupervised2'
os.makedirs(output_dir, exist_ok=True)

# 字体设置 - 标注字体25-30
LABEL_FONTSIZE = 26
TITLE_FONTSIZE = 28
TICK_FONTSIZE = 22
ANNOT_FONTSIZE = 18

sns.set_context("paper", font_scale=1.8)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.size'] = 14

# ==========================================
# 1. 核心计算函数
# ==========================================

def unsupervised_classification(df_all, h_vars, ele_map, n_clusters=5):
    """
    全数据标准化和评分，仅K-Means按高程分区执行。
    聚类以（高程区，cluster）唯一标识，避免不同分区的标签混淆。
    汇总所有分区的聚类均值，按全局分位数组建立统一分级阈值。
    返回分类结果（1-5，1=危害性最小，5=危害性最大）
    """
    # 1. 标准化：使用全部数据拟合一次
    scaler = StandardScaler()
    X = df_all[h_vars].fillna(0)
    X_scaled = pd.DataFrame(scaler.fit_transform(X), index=df_all.index, columns=h_vars)

    # 2. K-Means聚类：各高程分区独立拟合，输入使用全数据标准化结果
    df_temp = df_all.copy()
    df_temp['cluster'] = pd.Series(pd.NA, index=df_all.index, dtype='Int64')
    for ele_val in ele_map:
        mask = df_all['Altitude'] == ele_val
        if mask.sum() < max(10, n_clusters):
            continue
        kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        df_temp.loc[mask, 'cluster'] = kmeans.fit_predict(X_scaled.loc[mask])

    # 3. 计算每个类别的"危害性得分"（各因子Min-Max归一化值的均值）
    # Min-Max的最小值、最大值均来自全部数据；常量因子和缺失值计为0
    factor_min = df_all[h_vars].min()
    factor_range = (df_all[h_vars].max() - factor_min).replace(0, np.nan)
    norm_h = (df_all[h_vars] - factor_min) / factor_range
    norm_h = norm_h.fillna(0)
    df_temp['hazard_score'] = norm_h.mean(axis=1)

    # 在全部已聚类数据中，计算各（高程区，聚类）的平均得分，最多25个
    cluster_scores = df_temp.groupby(['Altitude', 'cluster'])['hazard_score'].mean().sort_values()

    # 4. 全局分位数分级：每个聚类均值权重相同，不按样本数加权。
    # 25个均值升序分为五组，每组5个；以相邻组边界得分中点为阈值。
    # 使用原始精度计算，边界值归入较高等级；不强制各分区拥有全部五级。
    sorted_scores = np.sort(cluster_scores.to_numpy())
    class_thresholds = np.array([], dtype=float)
    if len(sorted_scores) > 0:
        score_groups = np.array_split(sorted_scores, min(5, len(sorted_scores)))
        class_thresholds = np.array([
            (lower_group[-1] + upper_group[0]) / 2
            for lower_group, upper_group in zip(score_groups[:-1], score_groups[1:])
            if lower_group[-1] < upper_group[0]
        ], dtype=float)
        # 相邻组边界得分相同则合并该分界，避免相同得分被分为不同等级。
        if len(class_thresholds) < len(score_groups) - 1:
            print("分位数组边界存在相同得分，已合并对应分界。")
    cluster_classes = pd.Series(
        np.searchsorted(class_thresholds, cluster_scores.to_numpy(), side='right') + 1,
        index=cluster_scores.index, dtype=int)
    print(f"全局分位数分级使用的聚类数: {len(sorted_scores)}")
    print(f"危险性全局分位数分级阈值（相邻组边界中点）: {class_thresholds.tolist()}")
    rank_map = cluster_classes.to_dict()
    cluster_keys = pd.MultiIndex.from_frame(df_temp[['Altitude', 'cluster']])
    # 未聚类的少量样本沿用原脚本的0标记
    hazard_class = pd.Series(
        cluster_classes.reindex(cluster_keys).fillna(0).to_numpy(dtype=int),
        index=df_all.index)

    return hazard_class, cluster_scores, rank_map, df_temp['cluster']


def elevation_sensitivity_risk(df, h_vars, baseline_cluster_scores, grid_path=None,
                               return_h=False, return_elevation=False):
    """只重新计算分区方案的H，所有方案沿用输入数据的V和P。"""
    # 新增依赖：pip install rasterio pyshp shapely pyproj
    from pathlib import Path
    import rasterio
    from rasterio.mask import mask as raster_mask
    import shapefile
    from shapely.geometry import shape, mapping
    from shapely.ops import transform
    from pyproj import CRS, Transformer

    source_dir = Path(__file__).resolve().parent / '数据源'
    grid_path = Path(grid_path) if grid_path is not None else source_dir / '明长城选中1km.shp'
    dem_path = source_dir / 'dem.tif'
    if 'FID' not in df.columns or 'R' not in df.columns:
        raise KeyError('高程敏感性计算需要表格FID和原风险值R。')
    ids = pd.to_numeric(df['FID'], errors='raise')
    if ids.isna().any() or (ids != np.floor(ids)).any() or ids.duplicated().any():
        raise ValueError('表格FID必须是非空、唯一的整数。')
    grid_crs = CRS.from_wkt(grid_path.with_suffix('.prj').read_text())
    with shapefile.Reader(str(grid_path), encoding='utf-8') as reader:
        geometries = {
            feature_fid: shape(feature_shape.__geo_interface__)
            for feature_fid, feature_shape in enumerate(reader.iterShapes())
        }
    missing_ids = sorted(set(ids.astype(int)) - set(geometries))
    if missing_ids:
        raise ValueError(f'表格FID在网格中不存在: {missing_ids[:10]}')

    elevations = pd.Series(np.nan, index=df.index, dtype=float)
    with rasterio.open(dem_path) as dem:
        if dem.crs is None:
            raise ValueError('DEM缺少坐标系，无法安全提取高程。')
        project = Transformer.from_crs(grid_crs, dem.crs, always_xy=True)
        print(f'高程敏感性分析：网格CRS={grid_crs.to_string()}，DEM CRS={dem.crs}')
        for row_index, feature_fid in ids.items():
            geom = transform(project.transform, geometries[int(feature_fid)])
            try:
                pixels, _ = raster_mask(dem, [mapping(geom)], crop=True,
                                        indexes=1, filled=False, all_touched=False)
            except ValueError as exc:
                raise ValueError(f'FID={feature_fid}与DEM无交集。') from exc
            values = pixels.compressed().astype(float)
            values = values[np.isfinite(values)]
            if len(values) == 0:
                raise ValueError(f'FID={feature_fid}内没有有效DEM像元。')
            elevations.loc[row_index] = np.mean(values * dem.scales[0] + dem.offsets[0])
    print(f'DEM网格平均高程范围（m）: {elevations.min():.3f} 至 {elevations.max():.3f}')
    dem_500_zones = np.clip(np.floor(elevations / 500).astype(int) + 1, 1, 5)
    print(f'原Altitude与网格DEM均值五区一致: {(df["Altitude"] == dem_500_zones).sum()}/{len(df)}')

    # 沿用原H计算中的全局标准化和等权Min-Max得分。
    X = pd.DataFrame(StandardScaler().fit_transform(df[h_vars].fillna(0)),
                     index=df.index, columns=h_vars)
    factor_min = df[h_vars].min()
    factor_range = (df[h_vars].max() - factor_min).replace(0, np.nan)
    scores = ((df[h_vars] - factor_min) / factor_range).fillna(0).mean(axis=1)
    # 四种方案共用原始H的阈值，避免改变分区同时又改变分级标准。
    ordered = np.sort(baseline_cluster_scores.to_numpy())
    if not len(ordered):
        raise ValueError('原H没有有效聚类均值，无法建立共同分级阈值。')
    groups = np.array_split(ordered, min(5, len(ordered)))
    thresholds = np.array([(a[-1] + b[0]) / 2
                           for a, b in zip(groups[:-1], groups[1:]) if a[-1] < b[0]])
    print(f'四种分区方案共用原H阈值: {thresholds.tolist()}')
    result = pd.DataFrame(index=df.index)
    hazard_levels = pd.DataFrame(index=df.index)
    for width, column in [(250, 'R_250m'), (500, 'R_500m'),
                          (1000, 'R_1000m'), (None, 'R_no_elevation')]:
        hazard_column = column.replace('R_', 'H_', 1)
        if width == 500:
            # 原R已按Altitude的五区计算：0–500、500–1000、1000–1500、
            # 1500–2000和>2000 m。500 m基线应与原R完全相同。
            hazard_levels[hazard_column] = df['H']
            result[column] = df['R']
            print(f'{column}: 沿用原Altitude五区和R，{result[column].notna().sum()}个有效R值')
            continue
        # 250 m：0–2000 m每250 m一段，>=2000 m合为第九区。
        # 1000 m：0–1000、1000–2000、>=2000 m三个区。
        zones = (np.minimum(np.floor(elevations / width).astype(int), 2000 // width)
                 if width is not None else pd.Series(0, index=df.index))
        h = pd.Series(np.nan, index=df.index, dtype=float)
        for zone in sorted(zones.unique()):
            zone_index = zones.index[zones == zone]
            data = X.loc[zone_index]
            # 极小分区最多取其有效不同样本数，不把未聚类记录错误记为零风险。
            k = min(5, len(data.drop_duplicates()))
            labels = KMeans(n_clusters=k, random_state=42, n_init=10).fit_predict(data)
            cluster_means = scores.loc[zone_index].groupby(labels).mean()
            levels = pd.Series(np.searchsorted(thresholds, cluster_means, side='right') + 1,
                               index=cluster_means.index)
            h.loc[zone_index] = levels.reindex(labels).to_numpy()
            if k < 5:
                print(f'{column}: 分区{zone}仅可形成{k}个聚类（n={len(data)}）。')
        hazard_levels[hazard_column] = h
        result[column] = h * df['V'] / df['P'].replace(0, np.nan)
        print(f'{column}: {zones.nunique()}个分区，{result[column].notna().sum()}个有效R值')
    if return_elevation:
        return result, hazard_levels, elevations
    return (result, hazard_levels) if return_h else result


def plot_full_line_proportions(df, save_dir):
    """全线占比：H沿用现有等级，P用保护等级，V和R用Jenks自然断点。"""
    labels = ['Very Low', 'Low', 'Medium', 'High', 'Very High']
    colors = ['#38A800', '#8BD100', '#FFFF00', '#FF8000', '#FF0000']
    summaries = {}
    for col in ['H', 'V', 'P', 'R']:
        values = pd.to_numeric(df[col], errors='raise')
        valid = np.isfinite(values) & ((values >= 0) if col == 'V' else (values > 0))
        values = values.loc[valid]
        if values.empty:
            print(f'{col}: 无有效数据，跳过全线占比图。')
            continue

        if col == 'H':
            plot_labels, plot_colors = labels, colors
            if not values.isin(range(1, 6)).all():
                raise ValueError('H必须为1至5的整数等级。')
            counts = values.value_counts().reindex(range(1, 6), fill_value=0)
        elif col == 'P':
            # P=保护等级赋值/4*5，对应四种保护等级，不人为拆为五级。
            protection_scores = np.array([1.25, 2.5, 3.75, 5.0])
            if not np.isclose(values.to_numpy()[:, None], protection_scores).any(axis=1).all():
                raise ValueError('P与四种保护等级的预期赋值不一致。')
            plot_labels = ['No\ndesignation', 'Municipal/\ncounty', 'Provincial', 'National']
            plot_colors = ['#FF0000', '#FF8000', '#8BD100', '#38A800']
            categories = np.abs(values.to_numpy()[:, None] - protection_scores).argmin(axis=1)
            counts = pd.Series(categories).value_counts().reindex(range(4), fill_value=0)
        else:
            n_classes = min(5, values.nunique())
            if n_classes == 1:
                categories = np.zeros(len(values), dtype=int)
                plot_labels, plot_colors = ['Single class'], ['#FFFF00']
            else:
                breaks = jenkspy.jenks_breaks(values.to_numpy(), n_classes=n_classes)
                categories = np.searchsorted(breaks[1:-1], values.to_numpy(), side='left')
                print(f'{col}全线自然断点分级边界: {list(breaks)}')
                if n_classes == 5:
                    plot_labels, plot_colors = labels, colors
                else:
                    plot_labels = [f'Class {i + 1}' for i in range(n_classes)]
                    plot_colors = [colors[i] for i in np.rint(np.linspace(0, 4, n_classes)).astype(int)]
            counts = pd.Series(categories).value_counts().reindex(range(n_classes), fill_value=0)

        percentages = counts.to_numpy(dtype=float) / counts.sum() * 100
        fig, ax = plt.subplots(figsize=(12, 7))
        x = np.arange(len(plot_labels))
        bars = ax.bar(x, percentages, color=plot_colors, edgecolor='black', linewidth=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(plot_labels, fontsize=TICK_FONTSIZE)
        ax.set_xlabel('Protection Level' if col == 'P' else 'Level', fontsize=LABEL_FONTSIZE)
        ax.set_ylabel('Percentage (%)', fontsize=LABEL_FONTSIZE)
        ax.tick_params(axis='y', labelsize=TICK_FONTSIZE)
        ax.set_ylim(0, max(5, percentages.max() * 1.2))
        ax.set_axisbelow(True)
        ax.grid(axis='y', linestyle='--', alpha=0.4)
        offset = ax.get_ylim()[1] * 0.025
        for bar, percentage in zip(bars, percentages):
            ax.text(bar.get_x() + bar.get_width() / 2, percentage + offset,
                    f'{percentage:.2f}%', ha='center', va='bottom',
                    fontsize=28, fontweight='bold')
        fig.tight_layout()
        output_path = os.path.join(save_dir, f'{col}_Proportion_Plot.png')
        fig.savefig(output_path, bbox_inches='tight', dpi=300)
        plt.close(fig)
        summaries[col] = dict(zip(plot_labels, percentages))
        print(f'{col}全线占比图已保存: {output_path}（有效样本{len(values)}，排除{len(df) - len(values)}）')
    return summaries


def adaptive_discretization(series, col_name, n_bins=5):
    """自适应离散化：除了类别因子，全部使用自然断点法"""
    if col_name in ['CL', 'SA', 'Altitude']:
        return series

    n_unique = series.nunique()
    if n_unique < 2:
        return series

    value_counts = series.value_counts(normalize=True)
    if len(value_counts) > 0 and value_counts.iloc[0] > 0.95:
        return series

    actual_bins = min(n_bins, n_unique)

    try:
        # Jenks 自然断点
        breaks = jenkspy.jenks_breaks(series.dropna().values, n_classes=actual_bins)
        return pd.cut(series, bins=breaks, labels=False, include_lowest=True)
    except Exception:
        # 兜底：等频分箱（用rank避免重复边界）
        return pd.qcut(series.rank(method='first'), q=actual_bins, labels=False, duplicates='drop')


def factor_detector(df, feature_cols, target_col):
    """因子探测：计算解释力 q"""
    results = {}
    y = df[target_col]
    sst = y.var(ddof=0) * len(y)
    if sst == 0:
        return {col: 0.0 for col in feature_cols}

    for col in feature_cols:
        grouped = df.groupby(col)[target_col]
        ssw = sum(group.var(ddof=0) * len(group) for _, group in grouped)
        results[col] = max(0.0, 1.0 - (ssw / sst))
    return results


def get_interaction_type(q1, q2, q12):
    """判断交互类型"""
    if q12 < min(q1, q2):
        return "非线性削弱"
    if min(q1, q2) < q12 < max(q1, q2):
        return "单因子削弱"
    if q12 > (q1 + q2):
        return "非线性增强"
    if q12 > max(q1, q2):
        return "双因子增强"
    return "独立"


# ==========================================
# 2. 数据处理与计算
# ==========================================
file_path = 'RiskHEV8700.xlsx'
df_raw = pd.read_excel(file_path)

# 危害性因子
h_vars = ['EQ', 'LS', 'WE', 'WAE', 'FT', 'SA', 'NDVI', 'HRF', 'CL', 'RND', 'GI']
features = h_vars.copy()
ele_map = {1: 'ELER', 2: 'LER', 3: 'MER', 4: 'HER', 5: 'EHER'}

# 目标变量名（无监督分类结果）
target = 'H'

print("=" * 70)
print("按海拔分组进行无监督分类")
print("=" * 70)

# --- A. 按海拔进行无监督分类 ---
hazard_class_all, global_cluster_scores, global_rank_map, cluster_labels_all = unsupervised_classification(
    df_raw, h_vars, ele_map, n_clusters=5)
df_raw[target] = hazard_class_all

cluster_info = []         # 存储各海拔聚类信息
cluster_scores_all = {}   # 存储各海拔的聚类得分用于画图

for ele_val, ele_name in ele_map.items():
    mask = df_raw['Altitude'] == ele_val
    df_sub = df_raw[mask].copy()

    if len(df_sub) < 10:
        print(f"\n【{ele_name}】样本数不足({len(df_sub)}), 跳过")
        continue

    # 从全数据计算结果中提取本分区，用于汇总和绘图
    hazard_class = df_raw.loc[mask, target]
    cluster_scores = global_cluster_scores.xs(ele_val, level='Altitude').sort_values()
    rank_map = {k: global_rank_map[(ele_val, k)] for k in cluster_scores.index}

    # 记录聚类信息
    cluster_info.append({
        'Elevation': ele_name,
        'n': len(df_sub),
        'cluster_scores': cluster_scores.to_dict(),
        'cluster_classes': rank_map,
        'cluster_counts': cluster_labels_all.loc[mask].value_counts().to_dict(),
        'class_distribution': hazard_class.value_counts().sort_index().to_dict()
    })

    # 保留每个聚类的得分，避免多个聚类属于同一危害等级时相互覆盖
    cluster_scores_all[ele_name] = cluster_scores.to_dict()

    print(f"\n【{ele_name}】n={len(df_sub)}")
    print(f"  聚类危害性得分: {dict(cluster_scores.round(3))}")
    print(f"  聚类对应危害等级: {rank_map}")
    print(f"  分类分布: {hazard_class.value_counts().sort_index().to_dict()}")

# =========================
# 计算 V、P、R，连同危险性 H 共新增四列
# V = (Caozhizhi + Baohuzhuangtai / 6 * 5) / 2
# P = Baohudengjizhi / 4 * 5
# R = H * V / P
# =========================
required_cols = [target, 'Caozhizhi', 'Baohuzhuangtai', 'Baohudengjizhi']
missing = [c for c in required_cols if c not in df_raw.columns]
if missing:
    raise KeyError(f"缺少必要列：{missing}。请检查 RiskHEV8700.xlsx 的列名是否一致。")

material_vulnerability = pd.to_numeric(df_raw['Caozhizhi'], errors='raise')
preservation_status = pd.to_numeric(df_raw['Baohuzhuangtai'], errors='raise')
protection_level = pd.to_numeric(df_raw['Baohudengjizhi'], errors='raise')

df_raw['V'] = (material_vulnerability + preservation_status / 6 * 5) / 2
df_raw['P'] = protection_level / 4 * 5
df_raw['R'] = df_raw['H'] * df_raw['V'] / df_raw['P'].replace(0, np.nan)

# 三种替代风险公式；P_prime为P除以其理论最大值5，不新增辅助列。
alpha = 0.5  # 线性保护效应系数，可调整
beta = 1.0   # 指数保护效应系数，可调整
P_prime = df_raw['P'] / 5
df_raw['R2'] = df_raw['H'] * df_raw['V'] * (1 - alpha * P_prime)
df_raw['R3'] = df_raw['H'] * df_raw['V'] * np.exp(-beta * P_prime)
df_raw['R4'] = df_raw['H'] * df_raw['V']

# 以本表FID对应明长城选中1km.shp要素，从DEM取高程。
# 四种方案共用本表的V、P；500 m沿用原Altitude五区，使R_500m与R一致。
r_scales, h_scales, dem_mean = elevation_sensitivity_risk(
    df_raw, h_vars, global_cluster_scores, return_elevation=True)
scale_table = pd.concat([
    df_raw[['FID', 'Altitude', 'H', 'V', 'P', 'R']].reset_index(drop=True),
    dem_mean.rename('DEM_mean_m').reset_index(drop=True),
    h_scales.reset_index(drop=True), r_scales.reset_index(drop=True),
], axis=1)
scale_table.to_excel(f'{output_dir}/R_elevation_sensitivity.xlsx', index=False)
print('四种高程方案的H和R已保存至R_elevation_sensitivity.xlsx')

# 保存原始数据及 H、V、P、R 四列；缺失输入保持为缺失值
df_raw.to_excel(f'{output_dir}/data_with_hazard_class.xlsx', index=False)
print(f"\n带分类结果的数据已保存（已新增 H、V、P、R 四列）")

# 全线H、V、P、R统计图：不改变原始指数及既有分区图。
full_line_proportions = plot_full_line_proportions(df_raw, output_dir)

# --- B. 按海拔等级绘制因子排名图 ---
fig, axes = plt.subplots(2, 3, figsize=(26, 14))
axes = axes.flatten()

all_elevation_q = []
q_summary_list = []
interaction_list = []
risk_trends = []

for i, (ele_val, ele_name) in enumerate(ele_map.items()):
    df_sub = df_raw[df_raw['Altitude'] == ele_val].copy()
    if df_sub.empty or len(df_sub) < 10:
        continue

    # 离散化因子
    df_dis = df_sub.copy()
    for col in features:
        df_dis[col] = adaptive_discretization(df_sub[col], col)

    # 计算因子探测器 q 值（目标变量是无监督分类结果）
    current_q = factor_detector(df_dis, features, target)

    # 转换为 DataFrame 并排序
    rank_df = pd.DataFrame(list(current_q.items()), columns=['Feature', 'q']).sort_values('q', ascending=False)
    rank_df['Elevation'] = ele_name
    all_elevation_q.append(rank_df)

    # 绘图
    bars = sns.barplot(x='q', y='Feature', data=rank_df, ax=axes[i], palette='magma')
    axes[i].set_title(f'Elevation: {ele_name} (n={len(df_sub)})', fontsize=TITLE_FONTSIZE)
    axes[i].set_xlabel('Factor Determinant (q)', fontsize=LABEL_FONTSIZE)
    axes[i].set_ylabel('', fontsize=LABEL_FONTSIZE)
    axes[i].tick_params(axis='both', labelsize=TICK_FONTSIZE)

    # Reserve space inside each panel for all numeric labels.
    axes[i].set_xlim(0, max(float(rank_df['q'].max()) * 1.28, 0.1))
    for bar in bars.patches:
        width = bar.get_width()
        axes[i].annotate(f'{width:.3f}',
                         xy=(width, bar.get_y() + bar.get_height()/2),
                         xytext=(6, 0), textcoords='offset points',
                         ha='left', va='center', fontsize=ANNOT_FONTSIZE)

    # 存入 summary
    q_dict = current_q.copy()
    q_dict['Elevation_Range'] = ele_name
    q_summary_list.append(q_dict)

    # 风险探测（以 H 为“Mean Hazard Class”）
    for col in features:
        mean_risk = df_dis.groupby(col)[target].mean()
        for strata, val in mean_risk.items():
            risk_trends.append({'Elevation': ele_name, 'Factor': col, 'Strata': strata, 'Mean_Risk': val})

    # 交互探测（取 q 值前3的因子）
    top_factors = pd.Series(q_dict).drop('Elevation_Range').astype(float).sort_values(ascending=False).head(3).index.tolist()
    pairs = list(itertools.combinations(top_factors, 2))
    sst_sub = df_dis[target].var(ddof=0) * len(df_dis)

    for p1, p2 in pairs:
        df_dis['inter'] = df_dis[p1].astype(str) + "_" + df_dis[p2].astype(str)
        ssw_inter = sum(g.var(ddof=0) * len(g) for _, g in df_dis.groupby('inter')[target])
        q_inter = 1 - (ssw_inter / sst_sub if sst_sub != 0 else 0)
        interaction_list.append({
            'Elevation': ele_name, 'Factor_A': p1, 'Factor_B': p2,
            'q_A': q_dict[p1], 'q_B': q_dict[p2], 'q_Inter': q_inter,
            'Type': get_interaction_type(q_dict[p1], q_dict[p2], q_inter)
        })

# 移除多余的子图
for j in range(i + 1, len(axes)):
    fig.delaxes(axes[j])

fig.tight_layout(w_pad=3.0, h_pad=2.0)
fig.savefig(f'{output_dir}/factor_ranking_by_elevation.png', dpi=300)


# ==========================================
# 3. 结果汇总与图表输出
# ==========================================

# 保存 Excel 结果
q_summary_df = pd.DataFrame(q_summary_list).set_index('Elevation_Range')
q_summary_df.to_excel(f'{output_dir}/geodetector_q_summary.xlsx')

risk_df = pd.DataFrame(risk_trends)
risk_df.to_excel(f'{output_dir}/risk_detector.xlsx', index=False)

inter_df = pd.DataFrame(interaction_list)
inter_df.to_excel(f'{output_dir}/interaction_results.xlsx', index=False)

# 保存聚类信息
cluster_df = pd.DataFrame(cluster_info)
cluster_df.to_excel(f'{output_dir}/cluster_info.xlsx', index=False)

# ==========================================
# 4. 聚类信息可视化 (Cluster Info)
# ==========================================

# 4.1 各海拔聚类危害性得分对比图
fig, ax = plt.subplots(figsize=(14, 8))
ele_names = list(cluster_scores_all.keys())
x = np.arange(5)  # 每个分区的5个聚类，按平均危害性得分升序排列
width = 0.15
colors = sns.color_palette("Set2", len(ele_names))

for idx, ele_name in enumerate(ele_names):
    scores = cluster_scores_all[ele_name]
    values = list(scores.values())
    values += [np.nan] * (5 - len(values))
    bars = ax.bar(x + idx * width, values, width, label=ele_name, color=colors[idx])
    # 添加值标注
    for bar in bars:
        height = bar.get_height()
        if not np.isfinite(height):
            continue
        ax.text(bar.get_x() + bar.get_width()/2, height + 0.005,
                f'{height:.3f}', ha='center', va='bottom', fontsize=ANNOT_FONTSIZE, rotation=45)

ax.set_xlabel('Cluster Rank within Elevation', fontsize=LABEL_FONTSIZE)
ax.set_ylabel('Mean Hazard Score', fontsize=LABEL_FONTSIZE)
ax.set_title('Cluster Hazard Scores by Elevation', fontsize=TITLE_FONTSIZE)
ax.set_xticks(x + width * (len(ele_names) - 1) / 2)
ax.set_xticklabels(['Rank 1\n(Lowest)', 'Rank 2', 'Rank 3', 'Rank 4', 'Rank 5\n(Highest)'],
                   fontsize=TICK_FONTSIZE)
ax.tick_params(axis='y', labelsize=TICK_FONTSIZE)
ax.legend(title='Elevation', fontsize=TICK_FONTSIZE, title_fontsize=LABEL_FONTSIZE)
ax.set_ylim(0, max([max(v.values()) for v in cluster_scores_all.values()]) * 1.15)

plt.tight_layout()
plt.savefig(f'{output_dir}/cluster_hazard_scores.png', dpi=300)


# 4.2 各海拔类别分布堆叠柱状图
fig, ax = plt.subplots(figsize=(14, 8))
class_dist_data = {ele_name: [] for ele_name in ele_names}
for info in cluster_info:
    ele_name = info['Elevation']
    dist = info['class_distribution']
    class_dist_data[ele_name] = [dist.get(c, 0) for c in range(1, 6)]

# 转换为百分比
class_dist_pct = {}
for ele_name, counts in class_dist_data.items():
    total = sum(counts)
    class_dist_pct[ele_name] = [c / total * 100 if total > 0 else 0 for c in counts]

# 绘制堆叠柱状图
x = np.arange(len(ele_names))
bottom = np.zeros(len(ele_names))
colors = sns.color_palette("YlOrRd", 5)

for class_idx in range(5):
    values = [class_dist_pct[ele][class_idx] for ele in ele_names]
    bars = ax.bar(x, values, bottom=bottom, label=f'Class {class_idx + 1}', color=colors[class_idx])
    # 添加值标注
    for bar_idx, bar in enumerate(bars):
        height = bar.get_height()
        if height > 3:  # 只标注大于3%的
            ax.text(bar.get_x() + bar.get_width()/2, bottom[bar_idx] + height/2,
                    f'{height:.1f}%', ha='center', va='center', fontsize=ANNOT_FONTSIZE, color='black')
    bottom += values

ax.set_xlabel('Elevation', fontsize=LABEL_FONTSIZE)
ax.set_ylabel('Percentage (%)', fontsize=LABEL_FONTSIZE)
ax.set_title('Hazard Class Distribution by Elevation', fontsize=TITLE_FONTSIZE)
ax.set_xticks(x)
ax.set_xticklabels(ele_names, fontsize=TICK_FONTSIZE)
ax.tick_params(axis='y', labelsize=TICK_FONTSIZE)
ax.legend(title='Hazard Class', fontsize=TICK_FONTSIZE - 2, title_fontsize=LABEL_FONTSIZE - 2, loc='upper right')

plt.tight_layout()
plt.savefig(f'{output_dir}/hazard_class_distribution_stacked.png', dpi=300)


# ==========================================
# 5. 风险探测可视化 (Risk Detector)
# ==========================================

# 5.1 各因子在不同层级的平均风险热力图（按海拔分组）
for ele_name in ele_names:
    ele_risk = risk_df[risk_df['Elevation'] == ele_name]
    if ele_risk.empty:
        continue

    pivot_risk = ele_risk.pivot_table(index='Factor', columns='Strata', values='Mean_Risk', aggfunc='mean')

    fig, ax = plt.subplots(figsize=(14, 10))
    sns.heatmap(pivot_risk, annot=True, cmap='YlOrRd', fmt='.2f',
                annot_kws={'fontsize': ANNOT_FONTSIZE}, ax=ax,
                cbar_kws={'label': 'Mean Hazard Class'})
    ax.set_title(f'Risk Detector - {ele_name}', fontsize=TITLE_FONTSIZE)
    ax.set_xlabel('Strata (Discretized Level)', fontsize=LABEL_FONTSIZE)
    ax.set_ylabel('Factor', fontsize=LABEL_FONTSIZE)
    ax.tick_params(axis='both', labelsize=TICK_FONTSIZE)

    plt.tight_layout()
    plt.savefig(f'{output_dir}/risk_detector_{ele_name}.png', dpi=300)
    

# 5.2 主要因子的风险趋势折线图
top_factors_global = q_summary_df.mean().sort_values(ascending=False).head(6).index.tolist()

fig, axes = plt.subplots(2, 3, figsize=(20, 14))
axes = axes.flatten()

for idx, factor in enumerate(top_factors_global):
    ax = axes[idx]
    for ele_name in ele_names:
        factor_risk = risk_df[(risk_df['Elevation'] == ele_name) & (risk_df['Factor'] == factor)]
        if not factor_risk.empty:
            factor_risk_sorted = factor_risk.sort_values('Strata')
            ax.plot(factor_risk_sorted['Strata'], factor_risk_sorted['Mean_Risk'],
                    marker='o', label=ele_name, linewidth=2, markersize=8)
            # 添加值标注
            for xx, yy in zip(factor_risk_sorted['Strata'], factor_risk_sorted['Mean_Risk']):
                ax.annotate(f'{yy:.2f}', (xx, yy), textcoords="offset points",
                            xytext=(0, 8), ha='center', fontsize=ANNOT_FONTSIZE - 4)

    ax.set_title(f'Factor: {factor}', fontsize=TITLE_FONTSIZE)
    ax.set_xlabel('Strata', fontsize=LABEL_FONTSIZE)
    ax.set_ylabel('Mean Hazard Class', fontsize=LABEL_FONTSIZE)
    ax.tick_params(axis='both', labelsize=TICK_FONTSIZE)
    ax.legend(title='Elevation', fontsize=TICK_FONTSIZE - 4, title_fontsize=LABEL_FONTSIZE - 4)
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f'{output_dir}/risk_trend_top_factors.png', dpi=300)


# ==========================================
# 6. 各海拔分类分布图（单独柱状图）
# ==========================================
fig, axes = plt.subplots(1, 5, figsize=(25, 6))
for i, (ele_val, ele_name) in enumerate(ele_map.items()):
    df_sub = df_raw[df_raw['Altitude'] == ele_val]
    if len(df_sub) < 10:
        continue
    class_counts = df_sub[target].value_counts().sort_index()
    bars = axes[i].bar(class_counts.index, class_counts.values, color=sns.color_palette("YlOrRd", 5))
    axes[i].set_title(f'{ele_name}', fontsize=TITLE_FONTSIZE)
    axes[i].set_xlabel('Hazard Class', fontsize=LABEL_FONTSIZE)
    axes[i].set_ylabel('Count', fontsize=LABEL_FONTSIZE)
    axes[i].set_xticks([1, 2, 3, 4, 5])
    axes[i].tick_params(axis='both', labelsize=TICK_FONTSIZE)

    for bar in bars:
        height = bar.get_height()
        axes[i].text(bar.get_x() + bar.get_width()/2, height + max(class_counts) * 0.02,
                     f'{int(height)}', ha='center', va='bottom', fontsize=ANNOT_FONTSIZE)

plt.tight_layout()
plt.savefig(f'{output_dir}/hazard_class_distribution.png', dpi=300)

# 6.1 各海拔五个聚类的占比 + 原始聚类均值（双轴图）
# 按聚类均值升序排列；等级只决定横轴文字和柱色，同等级聚类不合并。
pct_labels = ['Very low', 'Low', 'Medium', 'High', 'Very high']
pct_colors = ['#38A800', '#8BD100', '#FFFF00', '#FF8000', '#FF0000']
fig_pct, axes_pct = plt.subplots(1, len(ele_map), figsize=(32, 10), squeeze=False)

for ax_pct, (ele_val, ele_name) in zip(axes_pct[0], ele_map.items()):
    ax_pct.set_title(f'Elevation: {ele_name}', fontsize=24, pad=20)
    pct_mask = (df_raw['Altitude'] == ele_val) & cluster_labels_all.notna()
    if not pct_mask.any():
        ax_pct.text(0.5, 0.5, 'No classified samples', transform=ax_pct.transAxes,
                    ha='center', va='center', fontsize=18)
        ax_pct.set_xticks([])
        continue
    pct_scores = global_cluster_scores.xs(ele_val, level='Altitude').sort_values()
    pct_counts = cluster_labels_all.loc[pct_mask].value_counts().reindex(pct_scores.index, fill_value=0)
    pct_total = pct_counts.sum()
    pct_values = pct_counts / pct_total * 100
    pct_levels = [global_rank_map[(ele_val, cluster)] for cluster in pct_scores.index]
    pct_x = np.arange(len(pct_scores))

    pct_bars = ax_pct.bar(pct_x, pct_values.values,
                          color=[pct_colors[level - 1] for level in pct_levels],
                          edgecolor='black', alpha=0.8)
    ax_pct.set_ylim(0, 80)
    ax_pct.set_ylabel('Percentage (%)', fontsize=20)
    ax_pct.set_xticks(pct_x)
    ax_pct.set_xticklabels([pct_labels[level - 1] for level in pct_levels],
                           fontsize=20, rotation=25)
    ax_pct.tick_params(axis='y', labelsize=20)
    for bar in pct_bars:
        inside = bar.get_height() >= 8
        label_y = bar.get_height() - 2 if inside else bar.get_height() + 1
        ax_pct.text(bar.get_x() + bar.get_width() / 2, label_y,
                    f'{bar.get_height():.0f}%', ha='center', va='top' if inside else 'bottom',
                    fontsize=18, fontweight='bold')

    ax_score = ax_pct.twinx()
    ax_score.plot(pct_x, pct_scores.values, color='black', marker='D',
                  markersize=10, linewidth=3, label='Mean Score')
    ax_score.set_ylabel('Mean Hazard Score', fontsize=20)
    ax_score.set_ylim(0, 0.4)
    ax_score.tick_params(axis='y', labelsize=20)
    for x_pos, score in zip(pct_x, pct_scores.values):
        if np.isfinite(score):
            ax_score.annotate(f'{score:.2f}', (x_pos, score),
                              textcoords='offset points', xytext=(0, 10),
                              ha='center', fontsize=18, fontweight='bold')
    if ele_val == next(iter(ele_map)):
        ax_score.legend(loc='upper left', fontsize=16)

fig_pct.tight_layout()
fig_pct.savefig(f'{output_dir}/hazard_class_distribution_pct.png', dpi=300)
plt.close(fig_pct)


# ==========================================
# 7. 因子q值汇总热力图
# ==========================================
fig, ax = plt.subplots(figsize=(16, 10))
q_heatmap = q_summary_df.T
sns.heatmap(q_heatmap, annot=True, cmap='YlOrRd', fmt='.3f',
            annot_kws={'fontsize': ANNOT_FONTSIZE}, ax=ax, cbar_kws={'label': 'q-value'})
ax.set_title('Factor Detector q-values by Elevation', fontsize=TITLE_FONTSIZE)
ax.set_xlabel('Elevation', fontsize=LABEL_FONTSIZE)
ax.set_ylabel('Factor', fontsize=LABEL_FONTSIZE)
ax.tick_params(axis='both', labelsize=TICK_FONTSIZE)

plt.tight_layout()
plt.savefig(f'{output_dir}/factor_q_heatmap.png', dpi=300)


# ==========================================
# 8. 交互探测热力图
# ==========================================
if len(inter_df) > 0:
    fig, ax = plt.subplots(figsize=(16, 12))
    pivot_inter = inter_df.pivot_table(index="Elevation", columns=["Factor_A", "Factor_B"], values="q_Inter")
    elevation_order = list(ele_map.values())
    pivot_inter = pivot_inter.reindex(elevation_order).dropna(how="all")
    sns.heatmap(pivot_inter, annot=True, cmap="YlOrRd", fmt=".3f",
                annot_kws={'fontsize': ANNOT_FONTSIZE}, ax=ax)
    ax.set_title('Interaction Detector Summary (q_Inter) - Unsupervised Classification', fontsize=TITLE_FONTSIZE)
    ax.set_xlabel('Factor Pairs', fontsize=LABEL_FONTSIZE)
    ax.set_ylabel('Elevation', fontsize=LABEL_FONTSIZE)
    ax.tick_params(axis='both', labelsize=TICK_FONTSIZE - 2)

    plt.tight_layout()
    plt.savefig(f'{output_dir}/interaction_heatmap_summary.png', dpi=300)
    

# ==========================================
# 9. 各海拔因子重要性雷达图
# ==========================================
from math import pi

fig, axes = plt.subplots(2, 3, figsize=(20, 14), subplot_kw=dict(polar=True))
axes = axes.flatten()

for idx, ele_name in enumerate(ele_names):
    if idx >= len(axes):
        break
    ax = axes[idx]
    if ele_name not in q_summary_df.index:
        continue

    q_vals = q_summary_df.loc[ele_name].values
    factors = q_summary_df.columns.tolist()
    N = len(factors)

    angles = [n / float(N) * 2 * pi for n in range(N)]
    angles += angles[:1]
    q_vals_plot = list(q_vals) + [q_vals[0]]

    ax.set_theta_offset(pi / 2)
    ax.set_theta_direction(-1)
    ax.set_xticks([n / float(N) * 2 * pi for n in range(N)])
    ax.set_xticklabels(factors, fontsize=TICK_FONTSIZE - 4)

    ax.plot(angles, q_vals_plot, 'o-', linewidth=2, color='steelblue')
    ax.fill(angles, q_vals_plot, alpha=0.25, color='steelblue')
    ax.set_title(f'{ele_name}', fontsize=TITLE_FONTSIZE, y=1.08)

for j in range(len(ele_names), len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout()
plt.savefig(f'{output_dir}/factor_radar_by_elevation.png', dpi=300)


# ==========================================
# 10. 样本量与平均风险汇总图
# ==========================================
fig, axes = plt.subplots(1, 2, figsize=(18, 7))

# 10.1 各海拔样本量
sample_counts = [info['n'] for info in cluster_info]
ele_labels = [info['Elevation'] for info in cluster_info]
bars1 = axes[0].bar(ele_labels, sample_counts, color=sns.color_palette("Blues", len(ele_labels)))
axes[0].set_xlabel('Elevation', fontsize=LABEL_FONTSIZE)
axes[0].set_ylabel('Sample Count', fontsize=LABEL_FONTSIZE)
axes[0].set_title('Sample Count by Elevation', fontsize=TITLE_FONTSIZE)
axes[0].tick_params(axis='both', labelsize=TICK_FONTSIZE)
for bar in bars1:
    height = bar.get_height()
    axes[0].text(bar.get_x() + bar.get_width()/2, height + max(sample_counts) * 0.02,
                 f'{int(height)}', ha='center', va='bottom', fontsize=ANNOT_FONTSIZE)

# 10.2 各海拔平均危害等级
mean_hazard = []
for ele_name in ele_labels:
    ele_val = [k for k, v in ele_map.items() if v == ele_name][0]
    mask = df_raw['Altitude'] == ele_val
    mean_hazard.append(df_raw.loc[mask, target].mean())

bars2 = axes[1].bar(ele_labels, mean_hazard, color=sns.color_palette("Oranges", len(ele_labels)))
axes[1].set_xlabel('Elevation', fontsize=LABEL_FONTSIZE)
axes[1].set_ylabel('Mean Hazard Class', fontsize=LABEL_FONTSIZE)
axes[1].set_title('Mean Hazard Class by Elevation', fontsize=TITLE_FONTSIZE)
axes[1].tick_params(axis='both', labelsize=TICK_FONTSIZE)
for bar in bars2:
    height = bar.get_height()
    axes[1].text(bar.get_x() + bar.get_width()/2, height + 0.05,
                 f'{height:.2f}', ha='center', va='bottom', fontsize=ANNOT_FONTSIZE)

plt.tight_layout()
plt.savefig(f'{output_dir}/sample_and_mean_hazard.png', dpi=300)


print(f"\n所有文件已生成并保存至: {output_dir}")


# ==========================================
# 11. 新增：q显著性检验与离散化敏感性分析（不修改以上既有计算）
# 结果单独写入 GeoDetector_Robustness 子目录。
# 单独补跑新增分析：python geodetector_robustness.py
# ==========================================
from geodetector_robustness import run_geodetector_robustness
geodetector_robustness_results = run_geodetector_robustness(
    df_raw, features, ele_map, target, adaptive_discretization,
    factor_detector, output_dir, permutations=9999, seed=42)
