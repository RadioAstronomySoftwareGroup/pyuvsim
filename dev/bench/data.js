window.BENCHMARK_DATA = {
  "lastUpdate": 1757359348081,
  "repoUrl": "https://github.com/RadioAstronomySoftwareGroup/pyuvsim",
  "entries": {
    "Benchmark": [
      {
        "commit": {
          "author": {
            "email": "burdorfmitchell@gmail.com",
            "name": "Mitchell Burdorf",
            "username": "burdorfmitchell"
          },
          "committer": {
            "email": "brynah@phys.washington.edu",
            "name": "Bryna Hazelton",
            "username": "bhazelton"
          },
          "distinct": true,
          "id": "52ace3f3baaa6b96a11f1e8997c63f16e90be2dc",
          "message": "Updated benchmark workflow logic to include pushes to main. Moved old benchmarks on gh-pages to old_bench.",
          "timestamp": "2025-09-08T12:16:10-07:00",
          "tree_id": "49fed0bf700103ac458d30c907ee5c2075973586",
          "url": "https://github.com/RadioAstronomySoftwareGroup/pyuvsim/commit/52ace3f3baaa6b96a11f1e8997c63f16e90be2dc"
        },
        "date": 1757359347032,
        "tool": "pytest",
        "benches": [
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.5_uvbeam]",
            "value": 0.04077254101447768,
            "unit": "iter/sec",
            "range": "stddev: 0.23725285906384885",
            "extra": "mean: 24.5263104805 sec\nrounds: 2"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_time_axis]",
            "value": 0.037218354641655925,
            "unit": "iter/sec",
            "range": "stddev: 0.05501867544632235",
            "extra": "mean: 26.868463413499995 sec\nrounds: 2"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.4_source_axis]",
            "value": 0.029866556101599873,
            "unit": "iter/sec",
            "range": "stddev: 0.06964966297928758",
            "extra": "mean: 33.482266806999974 sec\nrounds: 2"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_frequency_axis]",
            "value": 0.026426770831055525,
            "unit": "iter/sec",
            "range": "stddev: 0.26511228782934",
            "extra": "mean: 37.84041593250001 sec\nrounds: 2"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_baseline_number]",
            "value": 0.30389310916291007,
            "unit": "iter/sec",
            "range": "stddev: 0.10249789413472189",
            "extra": "mean: 3.2906307180000027 sec\nrounds: 2"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.7_multi_beam]",
            "value": 0.011814441097303542,
            "unit": "iter/sec",
            "range": "stddev: 0.5371065001847195",
            "extra": "mean: 84.64217577150002 sec\nrounds: 2"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.8_lunar]",
            "value": 0.30621648115546096,
            "unit": "iter/sec",
            "range": "stddev: 0.08196948119029276",
            "extra": "mean: 3.2656635469999955 sec\nrounds: 2"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.6_healpix]",
            "value": 1.7205427768081334,
            "unit": "iter/sec",
            "range": "stddev: 0.0022960563286886764",
            "extra": "mean: 581.211936999992 msec\nrounds: 2"
          }
        ]
      }
    ]
  }
}