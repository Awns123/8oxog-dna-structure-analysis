# Reproducibility

## 두 재현 수준

### Level 1 — 공개 pair table부터 핵심 결과 재현

목적은 DSSR 라이선스나 외부 네트워크 없이 최종 해석을 결정한 수치를 빠르게 검증하는 것이다.

```bash
python -m pip install -r requirements.txt
python run_analysis.py
python -m unittest discover -s tests -v
```

입력:

- `data/processed/reference_pairs_full_v1.csv`
- `data/processed/target_pairs_full_v1.csv`

출력:

- `results/generated/target_signed_six_distances.csv`
- `results/generated/matched_radial_direct_comparison.csv`
- `results/generated/direct_component_contributions.csv`
- `results/generated/leave_one_variable_out.csv`
- `results/generated/key_results.json`

`run_analysis.py`는 19개 검사를 통과하지 못하면 오류로 종료한다. 검사는 기준쌍 수, 표적 수,
잠근 핵심 수치 13개, 183D의 경험 순위, stretch 제외 시 두 위치의 방향 역전을 포함한다.

### Level 2 — RCSB 좌표부터 재현

필요 조건:

- Python 3.12
- Node.js 18 이상
- 인터넷 연결
- 별도로 받은 DSSR `v2.9.1-2026jul09`
- Python 의존성: `requirements.txt`

기준 전체 실행은 Windows 11용 DSSR 실행파일의 버전과 SHA-256을 모두 고정한다. 같은 버전명의
macOS/Linux 바이너리는 해시가 다르므로 이 정확 재현 경로가 의도적으로 거부한다. 다른 운영체제
실행은 별도 환경 민감도 분석으로 기록해야 한다.

```bash
python run_full_pipeline.py --dssr-path "/path/to/x3dna-dssr"
```

전체 경로는 다음 순서로 실행된다.

1. 잠근 22개 입력파일 다운로드와 SHA-256 대조
2. mmCIF 입력 QC
3. 동일 DSSR 옵션으로 production 실행
4. 독립 재실행 manifest 대조
5. `pairs[].bp_params` 추출과 염기쌍 방향 정규화
6. 정상 기준패널·표적 매핑 검증
7. signed-six 거리, radial/direct 비교
8. 변수, family, 품질, 통계 정의 민감도 분석
   - 결정품질 감사용 validation XML 21개는 2026-08-04 manifest의 크기·SHA-256과 강제 대조
9. 그림 재생성
10. Level 1 회귀검사

마지막 Level 1 검사는 방금 생성된 `data/pipeline_workspace/04_parsed_pairs`를 직접 입력으로
사용한다. 기존 `data/processed` snapshot만 다시 읽어 허위 PASS가 발생하지 않도록 연결했다.

원좌표와 DSSR JSON은 저장소에 커밋되지 않는다. 생성 위치는 각각
`data/pipeline_workspace/01_raw_mmcif/full_v1_2026-08-04/`와
`data/pipeline_workspace/03_dssr_json/`이다.
원 validation XML도 커밋하지 않으며, 잠근 계보는
`config/validation_manifest_2026-08-04.csv`에 보존한다. 공식 서버의 파일이 바뀌어 잠금 해시와
다르면 전체 실행은 오류로 중단한다.

원자료와 DSSR 출력의 무결성을 보호하기 위해 전체 실행은 기존 실행 폴더를 자동 덮어쓰지 않는다.
같은 경로에서 두 번째 전체 실행을 시작하지 말고 새 clone 또는 새 worktree에서 실행한다.

## 고정 실행환경

기준 실행은 다음과 같다.

- 실행일: 2026-08-04
- OS: Windows 11
- Python: 3.12.13
- NumPy: 2.3.5
- pandas: 3.0.1
- DSSR: v2.9.1-2026jul09
- DSSR executable SHA-256:
  `2fbf5dd32df8a66753486b00fdb2388e4188a81f9ed7251378b2fa889b82a7dc`
- DSSR options: `--more --json --nt-mapping=8OG:g`

현재 빠른 경로는 `requirements.txt`의 잠근 패키지 버전으로도 검증한다. 전체 기준 실행의 세부정보는
`config/runtime_manifest_full_v1.json`에 있다.

## 재현성의 의미

여기서 PASS는 공개 입력과 고정한 계산 규칙으로 같은 수치가 재생성된다는 뜻이다. 이는 표본의
대표성, 생물학적 인과관계 또는 모집단 일반화를 보증하지 않는다. 계산 재현성과 결론 타당성은
구분한다.

## 무결성 확인

공개 전에는 다음을 실행한다.

```bash
python scripts/preflight_publication.py
```

이 검사는 Python 문법, 금지한 바이너리·문서·압축파일, 개인 로컬 경로, 전형적인 비밀정보 패턴,
`SHA256SUMS.csv`의 파일 크기와 해시를 확인한다.
