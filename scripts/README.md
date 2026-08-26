# Script map

| 단계 | 스크립트 | 기능 |
|---|---|---|
| 입력 | `download_dssr_full_inputs_v1.js` | 잠근 PDB ID의 공식 mmCIF 다운로드·manifest 생성 |
| 입력 QC | `qc_dssr_full_inputs_v1.py` | 파일 해시·좌표 기본조건 검사 |
| DSSR | `run_dssr_full_v1.py` | 동일 옵션으로 모든 좌표 실행·JSON 보존 |
| 재실행 | `reconcile_dssr_rerun_manifest_v1.py` | 비수치 메타데이터와 수치 payload를 분리해 반복 실행 대조 |
| 추출 | `extract_compare_dssr_full_v1.py` | `bp_params` 추출, pair deduplication, 방향 정규화, 감사표 생성 |
| 주 분석 | `analyze_dssr_recalculated_v1.py` | 기준통계, 표적 거리, matched comparison |
| 지표 | `analyze_issue2_radial_percentage_v1.py` | radial 증가율과 direct distance의 의미 구분 |
| 변수 | `analyze_issue3_variable_scope_v1.py` | 3/6변수, 변수별 기여, leave-one-variable-out |
| 기준패널 | `analyze_issue5_reference_family_weighting_v1.py` | family weighting과 leave-one-family-out |
| pair state | `analyze_issue6_pair_state_v1.py` | 염기쌍 상태의 보조 구조 비교 |
| 품질 | `analyze_residual_reference_quality_v1.py` | 잔여 품질 필터 민감도 |
| 결정품질 | `crystal_quality_audit_v1.py` | 구조 메타데이터·해상도 조건 감사 및 2026-08-04 validation XML 잠금 해시 대조 |
| 통계 | `statistical_residual_audit_v1.py` | covariance, bootstrap, multiverse, stratification |
| 구조검색 | `search_rcsb_8og_independent_structures_v1.js` | RCSB 8OG 후보 검색 |
| 구조검색 | `screen_rcsb_protein_free_8og_pairs_v1.py` | protein-free 8OG pair 후속 선별 |
| 그림 | `build_dna_figures_v3.py` | 최종 그림 1–5 재생성 |
| 노트북 | `build_final_notebook.py` | 독자용 최종 노트북 생성 |
| 노트북 | `execute_notebook_stdlib.py` | 노트북 코드셀을 상단부터 실행하고 출력 저장 |
| 공개검사 | `preflight_publication.py` | 문법·경로·비밀정보·바이너리·해시 검사 |
| 무결성 | `build_sha256_manifest.py` | 공개 파일 전체의 크기·SHA-256 manifest 생성 |

일반 검토자는 루트의 `run_analysis.py`만 실행해도 최종 결론을 결정한 수치를 확인할 수 있다.
`run_full_pipeline.py`는 위 좌표 재현 스크립트를 순서대로 호출한다.
