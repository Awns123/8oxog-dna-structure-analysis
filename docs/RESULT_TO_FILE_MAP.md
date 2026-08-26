# Result-to-file map

| 주장 또는 수치 | 주 계산 코드 | 입력 | 생성 결과 | 독립 검사·감사 |
|---|---|---|---|---|
| 정상 기준 230쌍 = A:T 105 + G:C 125 | `run_analysis.py` | `data/processed/reference_pairs_full_v1.csv` | `results/generated/key_results.json` | `tests/test_analysis.py` |
| 주 표적 5개와 183D 대칭 감사행 | `extract_compare_dssr_full_v1.py` | DSSR JSON | `data/processed/target_pairs_full_v1.csv` | `data/quality/symmetry_audit_183D_full_v1.csv` |
| site 4 radial ΔD = +4.268721... | `run_analysis.py` | 기준·표적 pair CSV | `results/generated/matched_radial_direct_comparison.csv` | 잠근 회귀값 + unit test |
| site 9 radial ΔD = +2.307804... | `run_analysis.py` | 기준·표적 pair CSV | 같은 파일 | 잠근 회귀값 + unit test |
| direct D = 7.409744, 4.194300 | `run_analysis.py` | 같은 입력 | 같은 파일 | 잠근 회귀값 + unit test |
| stretch direct D² 비중 = 73.44%, 69.24% | `run_analysis.py` | 같은 입력 | `results/generated/direct_component_contributions.csv` | 잠근 회귀값 + unit test |
| stretch 제외 ΔD = −1.294563, −0.628268 | `run_analysis.py` | 같은 입력 | `results/generated/leave_one_variable_out.csv` | 부호 역전 assertion |
| 183D D = 1.370434, G:C 기준 24/125 이하 | `run_analysis.py` | G:C 기준·183D 표적 | `results/generated/target_signed_six_distances.csv` | rank 24 회귀검사 |
| 방향 변환 규칙 | `extract_compare_dssr_full_v1.py` | DSSR raw pair | `data/quality/orientation_audit_full_v1.csv` | `OrientationAuditTests` |
| 183D 대칭 두 행의 변환 후 동일성 | 같은 코드 | 183D assembly1 | `data/quality/symmetry_audit_183D_full_v1.csv` | `test_183d_symmetry_rows_collapse_after_orientation` |
| family weighting·LOFO | `analyze_issue5_reference_family_weighting_v1.py` | reference family map | `results/reference/sensitivity/` 관련 파일 | 저장된 matrix와 target metrics |
| covariance·bootstrap·multiverse | `statistical_residual_audit_v1.py` | pair table·family map | `results/reference/sensitivity/statistical_residual_audit_*` | 고정 seed와 원자료 표 |
| 추가 8OG 구조 검색 | JS search + Python screen | RCSB API/DSSR 결과 | `results/reference/search/` | search manifest와 manual classification |
| 그림 1–5 | `build_dna_figures_v3.py` | core·sensitivity 결과 | `results/figures/*.png` | 그림 3 입력 CSV 포함 |

## 읽는 순서

1. `README.md`에서 질문과 결론 범위를 확인한다.
2. `run_analysis.py`를 실행해 19/19 검사와 핵심 수치를 재생성한다.
3. 이 표에서 관심 결과의 row-level CSV와 감사표를 연다.
4. 변수 제외, 방향 변환, 183D 대칭 처리를 우선 검토한다.
5. 일반화 가능성은 `KNOWN_LIMITATIONS.md`와 확장 민감도 파일을 함께 본다.
