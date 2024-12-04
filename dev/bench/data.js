window.BENCHMARK_DATA = {
  "lastUpdate": 1733282494471,
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
      },
      {
        "commit": {
          "author": {
            "email": "90595017+burdorfmitchell@users.noreply.github.com",
            "name": "burdorfmitchell",
            "username": "burdorfmitchell"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "ec8d35a104844e98b963b9de86483fc30fe87f91",
          "message": "added dummy counter (#513)\n\n* added dummy counter\r\n\r\n* [pre-commit.ci] auto fixes from pre-commit.com hooks\r\n\r\nfor more information, see https://pre-commit.ci\r\n\r\n---------\r\n\r\nCo-authored-by: pre-commit-ci[bot] <66853113+pre-commit-ci[bot]@users.noreply.github.com>",
          "timestamp": "2024-12-02T20:10:39-05:00",
          "tree_id": "1e3256c5db3c34aac8d9400c743730095d33758e",
          "url": "https://github.com/RadioAstronomySoftwareGroup/pyuvsim/commit/ec8d35a104844e98b963b9de86483fc30fe87f91"
        },
        "date": 1733189473619,
        "tool": "pytest",
        "benches": [
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_uniform]",
            "value": 0.02171602026810698,
            "unit": "iter/sec",
            "range": "stddev: 0.462244264730384",
            "extra": "mean: 46.0489531532 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_uniform]",
            "value": 0.00150939424111593,
            "unit": "iter/sec",
            "range": "stddev: 0",
            "extra": "mean: 662.517434319 sec\nrounds: 1"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_mwa]",
            "value": 0.07791683012359271,
            "unit": "iter/sec",
            "range": "stddev: 0.385001886850954",
            "extra": "mean: 12.834197674799999 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_gauss]",
            "value": 0.001542399707051831,
            "unit": "iter/sec",
            "range": "stddev: 0",
            "extra": "mean: 648.3403721020001 sec\nrounds: 1"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_gauss]",
            "value": 0.022012491433470973,
            "unit": "iter/sec",
            "range": "stddev: 0.25992488717788736",
            "extra": "mean: 45.42875135339999 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_gauss]",
            "value": 0.6534141717500245,
            "unit": "iter/sec",
            "range": "stddev: 0.014260212053853307",
            "extra": "mean: 1.5304228822000028 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_uniform]",
            "value": 0.6459835923441396,
            "unit": "iter/sec",
            "range": "stddev: 0.02662337218104386",
            "extra": "mean: 1.5480269341999986 sec\nrounds: 5"
          }
        ]
      },
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
          "id": "aa866bc8049e119e1c6bf81e57f7330ae60a8ea7",
          "message": "cleaned up compare_post_benchmark.yaml a bit. now need to test running compare_post_benchmark using completion of another workflow (pull request and push)",
          "timestamp": "2024-12-03T20:11:16-05:00",
          "tree_id": "53af1b8e3593a378cfa0b2265de6079955fe1ca7",
          "url": "https://github.com/RadioAstronomySoftwareGroup/pyuvsim/commit/aa866bc8049e119e1c6bf81e57f7330ae60a8ea7"
        },
        "date": 1733275524198,
        "tool": "pytest",
        "benches": [
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_uniform]",
            "value": 0.02142095745533887,
            "unit": "iter/sec",
            "range": "stddev: 0.1329905607341969",
            "extra": "mean: 46.683254102199996 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_uniform]",
            "value": 0.0015404285061489643,
            "unit": "iter/sec",
            "range": "stddev: 0",
            "extra": "mean: 649.170017309 sec\nrounds: 1"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_mwa]",
            "value": 0.07949308373785116,
            "unit": "iter/sec",
            "range": "stddev: 0.026527205071192637",
            "extra": "mean: 12.579710749399993 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_gauss]",
            "value": 0.0015397566454737725,
            "unit": "iter/sec",
            "range": "stddev: 0",
            "extra": "mean: 649.4532775289999 sec\nrounds: 1"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_gauss]",
            "value": 0.022065716039523408,
            "unit": "iter/sec",
            "range": "stddev: 0.19610275311620018",
            "extra": "mean: 45.3191728838 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_gauss]",
            "value": 0.6570862466845493,
            "unit": "iter/sec",
            "range": "stddev: 0.024614256536471094",
            "extra": "mean: 1.521870234000005 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_uniform]",
            "value": 0.6442578746054148,
            "unit": "iter/sec",
            "range": "stddev: 0.01822887418511201",
            "extra": "mean: 1.5521734998000056 sec\nrounds: 5"
          }
        ]
      },
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
          "id": "7f8799475203541d8f7ba73827128caeee40c98c",
          "message": "updated approach to computing num_mismatched and fixed style",
          "timestamp": "2024-12-03T22:05:31-05:00",
          "tree_id": "0e24a00c8a229c20f53158799e57b10ac952a9ea",
          "url": "https://github.com/RadioAstronomySoftwareGroup/pyuvsim/commit/7f8799475203541d8f7ba73827128caeee40c98c"
        },
        "date": 1733282493947,
        "tool": "pytest",
        "benches": [
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_uniform]",
            "value": 0.02167177461891152,
            "unit": "iter/sec",
            "range": "stddev: 0.16758745066400074",
            "extra": "mean: 46.1429678734 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_uniform]",
            "value": 0.0014936097478566405,
            "unit": "iter/sec",
            "range": "stddev: 0",
            "extra": "mean: 669.518929851 sec\nrounds: 1"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_mwa]",
            "value": 0.08006802359578045,
            "unit": "iter/sec",
            "range": "stddev: 0.06577784859140037",
            "extra": "mean: 12.489380342999993 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_gauss]",
            "value": 0.0015076711557790513,
            "unit": "iter/sec",
            "range": "stddev: 0",
            "extra": "mean: 663.2746114209999 sec\nrounds: 1"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_gauss]",
            "value": 0.02179822819790662,
            "unit": "iter/sec",
            "range": "stddev: 0.3388608554315344",
            "extra": "mean: 45.87528816199999 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_gauss]",
            "value": 0.6334959157185419,
            "unit": "iter/sec",
            "range": "stddev: 0.020439600957893277",
            "extra": "mean: 1.5785421424000048 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_uniform]",
            "value": 0.6466466664181103,
            "unit": "iter/sec",
            "range": "stddev: 0.0379831977535706",
            "extra": "mean: 1.546439581200002 sec\nrounds: 5"
          }
        ]
      }
    ]
  }
}