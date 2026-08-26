# Code history and provenance

## 원칙

이 저장소는 두 층을 구분한다.

1. `historical/2026-05/`: 당시 폴더에서 확인된 파일을 내용 변경 없이 보존한 역사적 코드
2. 현재 루트와 `scripts/`: 2026년 8월의 재계산·민감도 분석을 공개 환경에서 실행할 수 있도록
   상대경로와 명령 인터페이스를 정리한 권장 코드

현재 정리본을 과거 파일로 소급하지 않으며, 과거 코드가 현재의 최종 결론을 그대로 산출했다고
주장하지 않는다.

## 2026년 5월 보존 코드

| 경로 | SHA-256 |
|---|---|
| `week1-core/01_download_mmcif_and_metadata.sh` | `375f9d2f2babcdeb04ac9cb2fdd61f56c0fba4f1cab182e46105e4f0f5e4b788` |
| `week1-core/01_fetch_and_screen.py` | `09b8021673086acad25c9f67c04454c96eaaf011986c89fda5e4f4a436a382bd` |
| `week1-expanded/01_download_expanded_mmcif_and_metadata.sh` | `f4be6a2f623f064c95d155b95ac82c6c9814608ed17ab7b1ee33e886052ff33d` |
| `week1-expanded/01_screen_expanded_mmcif.py` | `fa7a818d9ab033e66cc57bbfc6cd6af08e71da9c538721399e3a367a66ac63e1` |
| `week1-expanded/01_screen_expanded_mmcif_v3.py` | `47b2bdd592c4a78a203ca269a725bb091093af27b547f375835d97ae5c0b86fb` |
| `week2-analysis/02_parse_mmcif_parameters.py` | `ed6aef8bdfe2d3e93292a36c294a7de101e7b23d46c2499c9e2a22b74d493899` |
| `week2-analysis/02_parse_ndb_geometry_v1.py` | `e5f44df56f7fafab81d70f3d7cf7a3c1db5cf886ba545e051279c11d55c7ea5d` |
| `week2-analysis/03_week2_quality_stats_and_figures_v1.py` | `e1fddb1eb610190cc8712c434a82383ad58522ef52275a633423231e697f4a9b` |
| `week2-analysis/03_week2_finalize_scores_and_figures.py` | `2b1bb63a844a2bf3e93d4c41485455906db019f5ca67f87a92bda127449b934e` |

위 표의 해시는 원본에서 확인한 값이다. 전체 공개 폴더의 현재 해시는 `SHA256SUMS.csv`에서 다시
확인할 수 있다.

## 현재 권장 코드에서 적용한 공개성 개선

- 사용자별 절대경로를 저장소 기준 상대경로로 변경
- DSSR 실행파일 경로를 코드에 고정하지 않고 `--dssr-path`/`DSSR_EXE`로 주입
- 좌표 없는 빠른 재현 경로 `run_analysis.py` 추가
- 전체 좌표 경로 `run_full_pipeline.py` 추가
- 핵심 수치와 방향 규칙을 자동 검사하는 `tests/` 추가
- 실행 결과가 포함된 최종 노트북과 노트북 생성·실행 검사 스크립트 추가
- 개인정보·비밀정보·금지 바이너리·파일 해시를 검사하는 publication preflight 추가
- 결과 주장과 원자료·코드·검사 사이의 연결표 추가

## 제외한 항목

- DSSR 실행파일과 배포문서
- 원고·DOCX·PDF 및 문서 패키징 코드
- 개인 로컬 경로가 포함된 구 노트북
- 임시 의존성 폴더, 캐시, 중첩 ZIP
- 8-oxoG 연구와 직접 관련 없는 DNA 곡률·단백질 구조 코드
