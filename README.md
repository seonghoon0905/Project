# Fast Approximate String Matching via FM-Index

인간 염색체 스케일(대용량)의 유전체 데이터에서 허용 불일치(Mismatch)가 포함된 리드(Reads) 문자열의 매핑 위치를 탐색하기 위한 알고리즘 프로젝트

## 프로젝트 구조

이 저장소는 두 가지 개별 모델로 구성되어 알고리즘 최적화의 성능 한계를 대조 분석

### 1. `Benchmark/` (기준 모델)
이론적 원형에 입각하여 구현된 Baseline 모델

**알고리즘**: 순환 행렬(Rotation Matrix) 정렬, Naive LF-Mapping, 전체 역추적 탐색
 **특징**: $O(N^2)$ 공간 복잡도를 가지며, N > 40,000,000 규모의 데이터에서 Out-Of-Memory(OOM)가 발생하는 시스템적 한계를 지님

### 2. `Project/` (제안 최적화 모델)
거대 유전체 데이터를 메모리 초과 없이 $O(N)$ 공간 복잡도로 처리하기 위해 설계된 모델
*   **BWT 최적화**: Suffix Array + `PySAIS` 라이브러리의 SA-IS 알고리즘을 연동
*   **LF-Mapping 최적화**: 간격 $k=64$ 단위의 빈도 체크포인트(Occ Sample)를 적용
*   **FM-Index 최적화**: 접미사 배열의 1/32 샘플링(`sa_sample`) 구조 구현
*   **매핑 알고리즘**: 비둘기집 원리 기반 기법 도입

## 실행 환경 및 요구사항

이 프로젝트는 파이썬 3.10 이상의 환경에서 구동
`Project` 빌드를 위해 아래 외부 라이브러리 설치가 필수적

```bash
pip install -r requirements.txt
```

* `PySAIS`: C언어 기반의 O(N) SA-IS 접미사 배열 생성

## 실행 방법

데이터 구축 및 벤치마크 테스트는 다음 파이프라인으로 수행

1. **데이터 생성** (scripts 폴더)
   ```bash
   python scripts/generate_data.py
   ```
   UCSC 서버에서 Human Chr21 데이터를 다운로드하고, Illumina 장비의 기계적 에러율과 인간 유전체 SNP 확률(0.5%) 모델이 결합된 정밀 시뮬레이터(Empirical Simulator)를 통해 500만 개의 리드(Reads.txt)를 추출하여 `Inputs/` 디렉토리에 일괄 배치

2. **매핑 알고리즘 구동**
   ```bash
   cd Project
   python src/main.py ../Inputs/config.txt

   # 벤치마크 모델 실행 (주의, 명시적으로 실행 방법을 달아두었지만, 대용량 데이터가 램 위에 올라가므로 절대 구동 금지)
   cd ../Benchmark
   python src/main.py ../Inputs/config.txt
   ```

실행이 완료되면 루트의 `Outputs/` 디렉토리에 정렬 결과를 담은 표준 `result.sam` 파일이 도출
