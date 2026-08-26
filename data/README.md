# Data directory

## `processed/`

빠른 재현의 입력이다. DSSR `pairs[].bp_params`에서 추출하고 염기쌍 방향을 통일한 pair-level
CSV를 담는다.

- `reference_pairs_full_v1.csv`: 정상 기준 230쌍(A:T 105, G:C 125)
- `target_pairs_full_v1.csv`: 표적 원시 감사행 6개, 주 분석 역할 5개

## `quality/`

표적 매핑, 방향 변환, 183D 대칭 중복, 입력 좌표 QC를 검토할 수 있는 감사표다.
`validation_report_full_v1.md`는 2026-08-04 동일 DSSR 재계산 직후의 중간 검증 보고서이며,
당시의 세 변수 거리 서술을 포함한다. 최종 signed-six 결론은 루트 README와
`results/generated/`를 기준으로 한다.

## `pipeline_workspace/`

전체 좌표 재현 스크립트가 사용하는 고정 폴더 구조와 공개 가능한 설정·파생 결과의 snapshot이다.
원 mmCIF, DSSR JSON, 실행 로그는 `.gitignore`로 제외한다.

원자료 출처와 재배포 방침은 `../DATA_AND_SOFTWARE_NOTICE.md`를 참조한다.
