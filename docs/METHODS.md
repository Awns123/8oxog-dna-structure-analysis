# Methods

## 1. 구조와 분석단위

정상 기준패널은 미리 고정한 18개 B-DNA X-ray 구조다. DSSR이 canonical Watson–Crick으로
인식한 비변형 DNA pair 중 A:T와 G:C를 포함했다. 주 표적은 111D G:A와 178D 8-oxoG:A의
대응 site 4·9이며, 183D 8-oxoG:C는 별도 sanity case다.

분석행은 염기쌍 하나지만 독립 추론단위가 모두 염기쌍인 것은 아니다. 동일 구조의 여러 pair와
site 4·9의 군집성을 인정하며, 구조·family 단위 민감도 분석을 별도로 수행한다.

## 2. 좌표와 DSSR

공식 RCSB mmCIF를 그대로 보존하고 SHA-256을 기록했다. DSSR
`v2.9.1-2026jul09`를 모든 구조에 같은 옵션
`--more --json --nt-mapping=8OG:g`으로 적용했다. 원 DSSR JSON의
`pairs[].bp_params`를 사용하고 `pairs[].bp_simpleParams`는 사용하지 않았다.

여섯 변수는 다음 순서다.

\[
x=(\mathrm{shear},\mathrm{stretch},\mathrm{stagger},
\mathrm{buckle},\mathrm{propeller},\mathrm{opening})
\]

앞의 세 값 단위는 Å, 뒤의 세 값 단위는 degree다.

## 3. 방향 정규화

DSSR의 염기쌍 변수는 염기 순서에 따라 일부 부호가 달라진다. canonical pair와 표적을 각각
A→T, G→C, G→A, 8OG→A, 8OG→C 방향으로 통일했다. 역순 pair의 DSSR 구분자가
`M+N`이면 여섯 변수의 부호를 모두 바꾸고, `M−N`이면 shear와 buckle만 바꾼다.

원값, 변환값, 사용 규칙은 `data/quality/orientation_audit_full_v1.csv`에 행 단위로 남겼다.
자동검사는 두 규칙과 183D 대칭행의 변환 후 일치를 확인한다.

## 4. 정상 기준공간

A:T와 G:C 기준을 분리한다. pair group \(g\)의 변수 \(j\)에 대해 평균과 표본 표준편차를
다음처럼 계산한다.

\[
\mu_{gj}=\frac{1}{n_g}\sum_i x_{ij},\qquad
s_{gj}=\sqrt{\frac{1}{n_g-1}\sum_i(x_{ij}-\mu_{gj})^2}
\]

주 분석은 pair-equal 평균과 표준편차다. 동일 구조의 pair 의존성과 기준패널 선택 문제는
structure/family weighting, leave-one-out, 품질필터, robust/covariance 분석으로 점검한다.

## 5. radial distance와 direct distance

표적 \(x\)가 정상 기준 중심에서 얼마나 떨어졌는지는 signed-six radial distance로 나타낸다.

\[
D_{radial}(x)=\sqrt{\sum_{j=1}^{6}
\left(\frac{x_j-\mu_j}{s_j}\right)^2}
\]

111D와 178D의 대응 표적 사이 차이는 direct standardized distance로 별도 계산한다.

\[
D_{direct}(x^{178},x^{111})=\sqrt{\sum_{j=1}^{6}
\left(\frac{x^{178}_j-x^{111}_j}{s_j}\right)^2}
\]

두 값은 답하는 질문이 다르다. `D178 − D111`은 기준 중심에서의 반지름 차이이며,
`Ddirect`는 두 표적 사이의 직접 거리다. 반지름 증가율을 두 구조 전체의 차이나 생물학적 효과
크기로 해석하지 않는다.

## 6. 변수 기여와 민감도

direct distance 제곱에서 변수 \(j\)의 비중은 다음과 같다.

\[
100\times\frac{((x^{178}_j-x^{111}_j)/s_j)^2}{D_{direct}^2}
\]

또한 여섯 변수 중 하나씩 제외해 radial difference의 부호와 크기를 다시 계산한다. stretch 제외
시 두 위치에서 부호가 모두 역전된 것이 결론 축소의 핵심 근거다.

추가 분석은 absolute 3/6변수, covariance와 shrinkage Mahalanobis 거리, structure/family
weighting, leave-one-structure/family-out, 해상도·terminal·잔차 품질 필터, family bootstrap,
stratification과 multiverse를 포함한다.

## 7. 183D symmetry

183D의 표적 pair는 ASU에서 완결되지 않아 공식 biological assembly 1을 사용했다. 대칭으로
표현된 두 행을 원시 감사표에 남겼고, 방향 정규화 후 여섯 값이 같음을 확인한 뒤 대표행 하나만
통계에 사용했다.

## 8. 결론 게이트

분석은 다음을 구분한다.

- 지지: 두 표적 구조 모두 정상 기준에서 큰 pair-internal 이탈을 보인다.
- 지지: 이 matched case에서 178D의 추가 radial 차이는 stretch 중심이다.
- 기각: 큰 이탈 자체가 8-oxoG에 특이적이다.
- 기각: 178D가 변수 전반에서 일관되게 더 이탈한다.
- 판단 불가: 산화의 모집단 인과효과와 기능적 결과.
