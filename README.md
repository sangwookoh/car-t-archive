# CAR-T Publication Archive

ImmuneTherapy Lab 홈페이지의 CAR-T 문헌 아카이브. PubMed에서 CAR-T · CAAR-T · CAR-NK
관련 논문을 자동 수집하여 검색 가능한 정적 페이지로 공개한다.

- 수집·배포 전 과정 자동. 사람 개입 없음.
- 레코드는 **추가만** 되고 과거 레코드는 수정되지 않는다 (append-only).
- GitHub Actions + GitHub Pages. 공개 저장소 기준 비용 없음.

## 구조

```
scripts/fetch_pubmed.py      PubMed 수집기 (E-utilities, 표준 라이브러리만 사용)
.github/workflows/
  weekly.yml                 매주 토요일 자동 수집
  backfill.yml               과거 연도 수동 수집 (Actions 탭에서 실행)
data/
  manifest.json              연도 목록과 총계
  idx-<년도>.json            서지 메타데이터 (초록 제외) — 검색 대상
  abs-<년도>.json            초록 (펼칠 때만 로드)
  sum-<년도>.json            한국어 요약 (있는 논문만)
index.html                   렌더링 페이지. data/ 에서 전부 생성
```

페이지는 최근 2개 연도만 먼저 불러오고 나머지는 필요할 때 불러온다. 수만 건이 쌓여도
첫 화면 속도가 유지된다.

## 설치

1. 이 저장소를 **공개(public)** 로 만든다. 비공개 저장소는 Pages 발행에 유료 플랜이 필요하다.
2. Settings → Pages → Source 를 **Deploy from a branch**, 브랜치 `main` / 폴더 `/ (root)` 로 지정.
3. Settings → Actions → General → Workflow permissions 를 **Read and write permissions** 로 설정.
4. (선택) Settings → Secrets and variables → Actions 에 다음을 등록하면 수집 속도가 3배 이상 빨라진다.
   - `NCBI_API_KEY` — https://account.ncbi.nlm.nih.gov/settings/ 에서 무료 발급
   - `NCBI_EMAIL` — NCBI 권장 사항

## 과거 논문 백필

Actions 탭 → **Backfill by year** → Run workflow → 연도 구간 입력.

권장 순서 (한 번에 한 구간씩, 완료 후 다음 구간):

| 순서 | 구간 | 예상 건수 |
|---|---|---|
| 1 | 1987–2014 | 약 420 |
| 2 | 2015–2018 | 약 1,500 |
| 3 | 2019–2020 | 약 3,000 |
| 4 | 2021–2022 | 약 4,500 |
| 5 | 2023–2024 | 약 6,500 |
| 6 | 2025–2026 | 약 6,500 |

이미 저장된 PMID는 건너뛰므로 같은 구간을 다시 실행해도 중복이 생기지 않는다.
중단되었으면 같은 구간을 그대로 다시 실행하면 된다.

## 주간 수집

`weekly.yml` 이 매주 토요일 00:47 UTC(한국시간 09:47)에 실행된다.
수집 창은 **최근 45일**이다. PubMed는 논문을 뒤늦게 색인하는 경우가 있어서
7일 창으로는 그 논문들이 영구히 누락된다. PMID 중복 제거가 있으므로
창을 넓게 잡아도 비용은 늘지 않는다.

## 한국어 요약

`data/sum-<년도>.json` 에 `{"PMID": {"ko": "...", "week": "...", "src": "..."}}` 형식으로
저장한다. 요약이 있는 논문에는 목록에서 `한국어 요약` 배지가 붙고 본문이 함께 표시된다.
요약이 없는 논문도 검색·필터·링크는 모두 정상 동작한다.

이 파일은 수집기가 건드리지 않으므로, 주간 요약 생성 작업이 별도로 append 하면 된다.

## 홈페이지 임베드

Google Sites 편집기에서 삽입 → 포함(Embed) → URL 또는 삽입 코드에
Pages 주소(`https://<사용자명>.github.io/<저장소명>/`)를 넣는다. 한 번만 하면 되고,
이후 내용이 바뀌어도 Sites 쪽은 수정할 필요가 없다.

## 로컬 확인

```bash
python3 -m http.server 8000
# http://localhost:8000
```

`file://` 로 열면 브라우저 보안 정책 때문에 `data/` 를 읽지 못한다.

## 출처

서지 정보는 PubMed (NCBI E-utilities)에서 가져온다. 각 레코드는 PMID와 DOI 링크를 포함한다.
