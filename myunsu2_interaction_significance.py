"""Table 5 joint-q permutation tests / 表5联合q置换检验。

Standalone: python myunsu2_interaction_significance.py
Reads saved H from myunsu2.py; checks every factor q against existing Table 4.
All 55 pairs per region are tested before the original top-three pairs are shown.
Does not test the statistical significance of enhancement over an additive null.
"""
from pathlib import Path
import ast
import hashlib
import itertools
import json

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
FEATURES = ['EQ', 'LS', 'WE', 'WAE', 'FT', 'SA', 'NDVI', 'HRF', 'CL', 'RND', 'GI']
ELEVATIONS = {1: 'ELER', 2: 'LER', 3: 'MER', 4: 'HER', 5: 'EHER'}


def load_definitions():
    """Extract pure functions without executing the original GIS/output pipeline."""
    import jenkspy
    scope = dict(np=np, pd=pd, jenkspy=jenkspy)
    for filename, names in [
        ('myunsu2.py', {'adaptive_discretization'}),
        ('geodetector_robustness.py', {'q_stat', 'permutation_test', 'bh_adjust', 'interaction_type'}),
    ]:
        path = HERE / filename
        tree = ast.parse(path.read_text(encoding='utf-8-sig'))
        nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in names]
        if {n.name for n in nodes} != names:
            raise ValueError(f'Missing definitions / 缺少函数: {filename}')
        exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), 'exec'), scope)
    return scope


def create_figures(selected, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    for chinese in [False, True]:
        with plt.rc_context({'font.family': 'Microsoft YaHei' if chinese else 'Arial',
                             'font.size': 12, 'axes.unicode_minus': False,
                             'pdf.fonttype': 42, 'ps.fonttype': 42}):
            fig, axes = plt.subplots(1, 5, figsize=(17.8, 5.0), sharey=True)
            for ax, name in zip(axes, ELEVATIONS.values()):
                sub = selected.loc[selected.Elevation == name]
                colors = ['#44779B' if t == 'Bivariate enhancement' else '#D99A44' for t in sub.Type]
                bars = ax.bar(range(3), sub.q_inter, width=.68, color=colors,
                              edgecolor='#333333', linewidth=.6)
                for bar, row in zip(bars, sub.itertuples()):
                    ax.text(bar.get_x() + bar.get_width()/2, row.q_inter + .025,
                            f'{row.q_inter:.3f}{row.Significance}', ha='center', va='bottom', fontsize=12)
                ax.set_title(name, fontsize=16, pad=12)
                ax.set_xticks(range(3), [f'{a} × {b}' for a, b in zip(sub.Factor_A, sub.Factor_B)],
                              rotation=45, ha='right', rotation_mode='anchor')
                ax.set_ylim(0, 1)
                ax.set_yticks(np.arange(0, 1.01, .2))
                ax.set_axisbelow(True)
                ax.grid(axis='y', color='#dddddd', linewidth=.6)
                ax.spines[['top', 'right']].set_visible(False)
            axes[0].set_ylabel('联合解释力 $q_{inter}$' if chinese else 'Joint explanatory power, $q_{inter}$', fontsize=14)
            handles = [Patch(facecolor='#44779B', label='双因子增强' if chinese else 'Bivariate enhancement'),
                       Patch(facecolor='#D99A44', label='非线性增强' if chinese else 'Nonlinear enhancement')]
            fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.5, .088), ncol=2, frameon=False, fontsize=12)
            fig.text(.5, .025,
                     '*** BH校正后 p < 0.001（全部275组检验）；星号表示联合q的显著性，不表示增强类型的显著性。' if chinese else
                     '*** BH-adjusted p < 0.001 (275 tests). Stars refer to joint q, not the significance of enhancement type.',
                     ha='center', fontsize=10.5)
            fig.subplots_adjust(left=.065, right=.99, top=.90, bottom=.31, wspace=.20)
            stem = 'Fig15_interaction_significance_CN' if chinese else 'Fig15_interaction_significance_EN'
            fig.savefig(out / f'{stem}.png', dpi=400, facecolor='white')
            fig.savefig(out / f'{stem}.pdf', facecolor='white')
            plt.close(fig)


def run_interaction_significance(df=None, output_dir=None, permutations=9999, seed=20261016):
    base = Path(output_dir) if output_dir else HERE / 'Results_Great/Hazard_Unsupervised2'
    out = base / 'Table5_Significance_20261016'
    out.mkdir(parents=True, exist_ok=True)
    funcs = load_definitions()
    if df is None:
        df = pd.read_excel(base / 'data_with_hazard_class.xlsx')
        raw = pd.read_excel(HERE / 'RiskHEV8700.xlsx')
        cols = ['FID', 'Altitude'] + FEATURES
        raw_matches = (raw[cols].shape == df[cols].shape and
                       np.array_equal(raw[cols].to_numpy(), df[cols].to_numpy(), equal_nan=True))
        if not raw_matches:
            print(f'NOTICE / 注意: current raw n={len(raw)}, saved Table4 baseline n={len(df)}; '
                  'using the Table4 baseline / 使用表4基线，不混用版本', flush=True)
    if not np.isfinite(df[['Altitude', 'H'] + FEATURES].to_numpy(float)).all():
        raise ValueError('Nonfinite inputs / 存在缺失或非有限输入')
    reference = pd.read_csv(base / 'GeoDetector_Robustness/q_significance_details.csv')
    expected_fingerprint = json.loads((base / 'GeoDetector_Robustness/analysis_settings.json').read_text(encoding='utf-8'))['input_fingerprint']
    fingerprint = hashlib.sha256(pd.util.hash_pandas_object(df[['Altitude', 'H'] + FEATURES], index=True).values.tobytes()).hexdigest()
    # Excel round-trip can change float bits and therefore pandas' exact hash.
    # The substantive gate below reproduces all 55 Table4 q values AND n,
    # then all 15 original interaction q values to absolute tolerance 1e-12.
    fingerprint_matches = fingerprint == expected_fingerprint
    if not fingerprint_matches:
        print('Exact binary fingerprint differs; checking all Table4/Table5 statistics / '
              '二进制指纹不同，将逐项核对全部表4与表5统计量', flush=True)
    old = pd.read_excel(base / 'interaction_results.xlsx')
    records = []
    for elevation, name in ELEVATIONS.items():
        sub = df.loc[df.Altitude == elevation].copy()
        y = sub.H.to_numpy(float)
        codes, q = {}, {}
        for factor in FEATURES:
            labels = funcs['adaptive_discretization'](sub[factor], factor)
            if labels.isna().any():
                raise ValueError(f'Missing discretized labels / 离散化标签缺失: {name}/{factor}')
            codes[factor] = pd.factorize(labels, sort=True)[0]
            q[factor] = funcs['q_stat'](y, codes[factor])
            ref_q = reference.loc[(reference.Elevation == name) & (reference.Factor == factor), 'q'].item()
            ref_n = reference.loc[(reference.Elevation == name) & (reference.Factor == factor), 'n'].item()
            if len(y) != ref_n:
                raise AssertionError(f'Table4 sample size mismatch / 表4样本量不同: {name}')
            if not np.isclose(q[factor], ref_q, rtol=0, atol=1e-12):
                raise AssertionError(f'Table 4 q mismatch / 表4不一致: {name}/{factor}')
        top = pd.Series(q).sort_values(ascending=False, kind='stable').head(3).index.tolist()
        selected_keys = {frozenset(p) for p in itertools.combinations(top, 2)}
        for pair_index, (a, b) in enumerate(itertools.combinations(FEATURES, 2)):
            joint_codes = pd.factorize(pd.MultiIndex.from_arrays([codes[a], codes[b]]), sort=True)[0]
            test_seed = seed + elevation * 1000 + pair_index
            q_inter, p, exceed = funcs['permutation_test'](y, joint_codes, permutations, test_seed)
            counts = np.bincount(joint_codes)
            # Independently calculate pooled within-stratum SS, matching myunsu2.py.
            within = sum(np.var(y[joint_codes == c], ddof=0) * n for c, n in enumerate(counts))
            direct = 1 - within / np.sum((y-y.mean())**2)
            if not np.isclose(q_inter, direct, rtol=0, atol=1e-12):
                raise AssertionError('Independent q check failed / 独立q核对失败')
            records.append(dict(Elevation=name, Factor_A=a, Factor_B=b, n=len(y),
                                q_A=q[a], q_B=q[b], q_inter=q_inter,
                                Type=funcs['interaction_type'](q[a], q[b], q_inter),
                                p_raw=p, Exceedances=exceed, Permutations=permutations, Seed=test_seed,
                                Joint_strata=len(counts), Min_joint_stratum_n=int(counts.min()),
                                Singleton_strata=int((counts==1).sum()),
                                Selected_for_Table5=frozenset((a,b)) in selected_keys))
        print(f'{name}: n={len(y)}, 55 pairs completed / 55组完成', flush=True)
    full = pd.DataFrame(records)
    if len(full) != 275:
        raise AssertionError('Expected 275 tests / 应为275项检验')
    full['p_BH'] = funcs['bh_adjust'](full.p_raw)
    full['Significance'] = full.p_BH.map(lambda p: '***' if p < .001 else '**' if p < .01 else '*' if p < .05 else 'ns')
    selected_rows = []
    for row in old.itertuples():
        hits = full[(full.Elevation == row.Elevation) &
                    (((full.Factor_A == row.Factor_A) & (full.Factor_B == row.Factor_B)) |
                     ((full.Factor_A == row.Factor_B) & (full.Factor_B == row.Factor_A)))]
        record = hits.iloc[0].to_dict()
        if not record['Selected_for_Table5']:
            raise AssertionError('Original selected pair mismatch / 原表筛选组合不一致')
        qa = record['q_A'] if record['Factor_A'] == row.Factor_A else record['q_B']
        qb = record['q_B'] if record['Factor_B'] == row.Factor_B else record['q_A']
        if not np.allclose([qa,qb,record['q_inter']], [row.q_A,row.q_B,row.q_Inter], atol=1e-12, rtol=0):
            raise AssertionError('Original Table 5 mismatch / 原表5的q值不一致')
        record.update(Factor_A=row.Factor_A, Factor_B=row.Factor_B, q_A=qa, q_B=qb)
        selected_rows.append(record)
    selected = pd.DataFrame(selected_rows)
    # Secondary audit only; manuscript uses the broader 275-test correction.
    selected['p_BH_15_audit_only'] = funcs['bh_adjust'](selected.p_raw)
    full.to_csv(out / 'all_275_joint_q_tests.csv', index=False, encoding='utf-8-sig')
    (out / 'all_275_joint_q_tests.json').write_text(full.to_json(orient='records', force_ascii=False, indent=2), encoding='utf-8')
    selected.to_csv(out / 'Table5_significance_details.csv', index=False, encoding='utf-8-sig')
    metadata = dict(permutations=permutations, seed=seed, target='H', n=len(df),
                    p_formula='(1 + count(q_perm >= q_observed))/(B + 1)',
                    correction='Benjamini-Hochberg across all 275 elevation-pair tests',
                    null='Exchangeable H labels within each elevation region; joint strata fixed',
                    selection='Original top-three-factor pairs shown; all 55 pairs per region tested before BH',
                    fingerprint=fingerprint, table4_recorded_fingerprint=expected_fingerprint,
                    binary_fingerprint_matches=fingerprint_matches,
                    verification='All 55 Table4 q and regional n, and all 15 Table5 qA/qB/q_inter reproduced at absolute tolerance 1e-12',
                    original_table4_and_table5_reproduced=True,
                    interpretation_en='Joint association test, not an additive/no-interaction test. Unrestricted permutations do not preserve spatial autocorrelation. H is constructed from these same factors; significance is not independent validation or causal evidence.',
                    interpretation_zh='检验联合空间关联，不检验相对于可加模型的交互增益。自由置换不保留空间自相关。H由这些因子构建，显著性不代表独立验证或因果证据。')
    (out / 'analysis_settings.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
    (out / 'workbook_data.json').write_text(selected.to_json(orient='records', force_ascii=False, indent=2), encoding='utf-8')
    create_figures(selected, out)
    print(selected[['Elevation','Factor_A','Factor_B','q_inter','Type','p_raw','p_BH','Significance']].to_string(index=False), flush=True)
    return selected


if __name__ == '__main__':
    run_interaction_significance()
