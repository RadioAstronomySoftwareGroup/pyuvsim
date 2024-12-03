window.BENCHMARK_DATA = {
  "lastUpdate": 1733185648329,
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
            "email": "burdorfmitchell@gmail.com",
            "name": "Mitchell Burdorf",
            "username": "burdorfmitchell"
          },
          "distinct": true,
          "id": "04c02f12d526d4f890fca4efe50e8011c0509068",
          "message": "commented out the 0 tolerance sim comparison check",
          "timestamp": "2024-12-02T19:13:50-05:00",
          "tree_id": "0f2f0292d8153493d85c752c97b80778c8440b3c",
          "url": "https://github.com/RadioAstronomySoftwareGroup/pyuvsim/commit/04c02f12d526d4f890fca4efe50e8011c0509068"
        },
        "date": 1733185648022,
        "tool": "pytest",
        "benches": [
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_uniform]",
            "value": 0.02202461193566478,
            "unit": "iter/sec",
            "range": "stddev: 0.238976935333406",
            "extra": "mean: 45.40375117260001 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_uniform]",
            "value": 0.0014947421494588043,
            "unit": "iter/sec",
            "range": "stddev: 0",
            "extra": "mean: 669.011709051 sec\nrounds: 1"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_mwa]",
            "value": 0.07976508245435068,
            "unit": "iter/sec",
            "range": "stddev: 0.016607325936829204",
            "extra": "mean: 12.536813969600008 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_gauss]",
            "value": 0.0015054060398427584,
            "unit": "iter/sec",
            "range": "stddev: 0",
            "extra": "mean: 664.272610534 sec\nrounds: 1"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_gauss]",
            "value": 0.020711289575964454,
            "unit": "iter/sec",
            "range": "stddev: 0.26298854455287485",
            "extra": "mean: 48.2828457558 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_gauss]",
            "value": 0.639482805542522,
            "unit": "iter/sec",
            "range": "stddev: 0.02080966684181971",
            "extra": "mean: 1.5637637029999951 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_uniform]",
            "value": 0.6567570612133533,
            "unit": "iter/sec",
            "range": "stddev: 0.02076744966992017",
            "extra": "mean: 1.5226330390000045 sec\nrounds: 5"
          }
        ]
      }
    ]
  }
}