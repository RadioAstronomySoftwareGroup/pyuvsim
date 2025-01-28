window.BENCHMARK_DATA = {
  "lastUpdate": 1738090168328,
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
        "date": 1733287925594,
        "tool": "pytest",
        "benches": [
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_uniform]",
            "value": 0.02142697109842286,
            "unit": "iter/sec",
            "range": "stddev: 0.14593284455581707",
            "extra": "mean: 46.6701520904 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_uniform]",
            "value": 0.0015220867544172745,
            "unit": "iter/sec",
            "range": "stddev: 0",
            "extra": "mean: 656.9927746220001 sec\nrounds: 1"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_mwa]",
            "value": 0.07833884957763176,
            "unit": "iter/sec",
            "range": "stddev: 0.015178368203211138",
            "extra": "mean: 12.765058529600003 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_gauss]",
            "value": 0.0015549797995768244,
            "unit": "iter/sec",
            "range": "stddev: 0",
            "extra": "mean: 643.095170929 sec\nrounds: 1"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_gauss]",
            "value": 0.021892307249642338,
            "unit": "iter/sec",
            "range": "stddev: 0.12154505638365004",
            "extra": "mean: 45.67814568820001 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_gauss]",
            "value": 0.6517244653205481,
            "unit": "iter/sec",
            "range": "stddev: 0.018707604381485114",
            "extra": "mean: 1.534390763600004 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_uniform]",
            "value": 0.647318015046799,
            "unit": "iter/sec",
            "range": "stddev: 0.022326276421076253",
            "extra": "mean: 1.5448357326000006 sec\nrounds: 5"
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
          "id": "738b8135d90d80410b7cb3a35cdb5794c83d14d7",
          "message": "not sure why linking to tests isn't working -- swapping back",
          "timestamp": "2024-12-04T00:58:58-05:00",
          "tree_id": "59755879d77ecdb19ca410479b121808ada138d6",
          "url": "https://github.com/RadioAstronomySoftwareGroup/pyuvsim/commit/738b8135d90d80410b7cb3a35cdb5794c83d14d7"
        },
        "date": 1733292845159,
        "tool": "pytest",
        "benches": [
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_uniform]",
            "value": 0.021365502515988285,
            "unit": "iter/sec",
            "range": "stddev: 0.20060450450924827",
            "extra": "mean: 46.80442218720002 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_uniform]",
            "value": 0.0015250150998164201,
            "unit": "iter/sec",
            "range": "stddev: 0",
            "extra": "mean: 655.731212183 sec\nrounds: 1"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_mwa]",
            "value": 0.07882749923829759,
            "unit": "iter/sec",
            "range": "stddev: 0.037783191118628275",
            "extra": "mean: 12.685928256800002 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_gauss]",
            "value": 0.0015451873096959797,
            "unit": "iter/sec",
            "range": "stddev: 0",
            "extra": "mean: 647.170730516 sec\nrounds: 1"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_gauss]",
            "value": 0.021625153111335702,
            "unit": "iter/sec",
            "range": "stddev: 0.056524716593715794",
            "extra": "mean: 46.24244715639999 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_gauss]",
            "value": 0.6463433692794459,
            "unit": "iter/sec",
            "range": "stddev: 0.017735864624275084",
            "extra": "mean: 1.5471652491999976 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_uniform]",
            "value": 0.6355515944833545,
            "unit": "iter/sec",
            "range": "stddev: 0.022321951556262667",
            "extra": "mean: 1.573436379799989 sec\nrounds: 5"
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
          "id": "ac59a8c0ec8da140f29f802ff64d14959512bb6b",
          "message": "edited README / environment.yaml to discuss/require requests and pytest-benchmark",
          "timestamp": "2024-12-04T02:28:02-05:00",
          "tree_id": "7a2df074d20b31c4313c34747a1ca4d63a73fd46",
          "url": "https://github.com/RadioAstronomySoftwareGroup/pyuvsim/commit/ac59a8c0ec8da140f29f802ff64d14959512bb6b"
        },
        "date": 1733298202350,
        "tool": "pytest",
        "benches": [
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_uniform]",
            "value": 0.021586843335568772,
            "unit": "iter/sec",
            "range": "stddev: 0.10623921268586099",
            "extra": "mean: 46.32451278099999 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_uniform]",
            "value": 0.001502531233523922,
            "unit": "iter/sec",
            "range": "stddev: 0",
            "extra": "mean: 665.543569204 sec\nrounds: 1"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_mwa]",
            "value": 0.0788661480047967,
            "unit": "iter/sec",
            "range": "stddev: 0.08193448608124185",
            "extra": "mean: 12.679711451599985 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_gauss]",
            "value": 0.0015034947539261823,
            "unit": "iter/sec",
            "range": "stddev: 0",
            "extra": "mean: 665.1170530450001 sec\nrounds: 1"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_gauss]",
            "value": 0.021299083097370435,
            "unit": "iter/sec",
            "range": "stddev: 0.15467194111142002",
            "extra": "mean: 46.95037788380004 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_gauss]",
            "value": 0.6367996289522346,
            "unit": "iter/sec",
            "range": "stddev: 0.021935985957498887",
            "extra": "mean: 1.5703526738000164 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_uniform]",
            "value": 0.6527405584277286,
            "unit": "iter/sec",
            "range": "stddev: 0.016005613343320128",
            "extra": "mean: 1.5320022436000045 sec\nrounds: 5"
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
          "id": "7725475b09bf36af03c2cb42b01cfcf53b229a93",
          "message": "edited README",
          "timestamp": "2024-12-04T21:58:14-05:00",
          "tree_id": "b617510d7974d5ed106cdd941463bd1a9774ca25",
          "url": "https://github.com/RadioAstronomySoftwareGroup/pyuvsim/commit/7725475b09bf36af03c2cb42b01cfcf53b229a93"
        },
        "date": 1733369717603,
        "tool": "pytest",
        "benches": [
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_uniform]",
            "value": 0.02165526009487669,
            "unit": "iter/sec",
            "range": "stddev: 0.0539076257273444",
            "extra": "mean: 46.17815697519999 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_uniform]",
            "value": 0.0015151563553460405,
            "unit": "iter/sec",
            "range": "stddev: 0",
            "extra": "mean: 659.997891618 sec\nrounds: 1"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_mwa]",
            "value": 0.07807296025851874,
            "unit": "iter/sec",
            "range": "stddev: 0.2584830748577664",
            "extra": "mean: 12.808531874399979 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_gauss]",
            "value": 0.0015267124956032854,
            "unit": "iter/sec",
            "range": "stddev: 0",
            "extra": "mean: 655.0021715810001 sec\nrounds: 1"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_gauss]",
            "value": 0.022157317457406144,
            "unit": "iter/sec",
            "range": "stddev: 0.12244274192488064",
            "extra": "mean: 45.131817149000014 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_gauss]",
            "value": 0.6349938694288279,
            "unit": "iter/sec",
            "range": "stddev: 0.024282740630429935",
            "extra": "mean: 1.5748183536 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_uniform]",
            "value": 0.6420379072111073,
            "unit": "iter/sec",
            "range": "stddev: 0.01854857003625209",
            "extra": "mean: 1.557540432999997 sec\nrounds: 5"
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
          "id": "b0e27c8f77e64592ef088e390cbaf06bec5667ff",
          "message": "edited README",
          "timestamp": "2024-12-04T22:26:26-05:00",
          "tree_id": "576b294565d7188b2ca7b522df78d30822bbbae4",
          "url": "https://github.com/RadioAstronomySoftwareGroup/pyuvsim/commit/b0e27c8f77e64592ef088e390cbaf06bec5667ff"
        },
        "date": 1733370661927,
        "tool": "pytest",
        "benches": [
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_uniform]",
            "value": 0.02170909806465009,
            "unit": "iter/sec",
            "range": "stddev: 0.11990961341351769",
            "extra": "mean: 46.0636364082 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_uniform]",
            "value": 0.0015000437311739063,
            "unit": "iter/sec",
            "range": "stddev: 0",
            "extra": "mean: 666.6472311560001 sec\nrounds: 1"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_mwa]",
            "value": 0.07937581658848629,
            "unit": "iter/sec",
            "range": "stddev: 0.03184332843781786",
            "extra": "mean: 12.598295589000003 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_gauss]",
            "value": 0.0015303675846535126,
            "unit": "iter/sec",
            "range": "stddev: 0",
            "extra": "mean: 653.4377818950001 sec\nrounds: 1"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_gauss]",
            "value": 0.02166064879342106,
            "unit": "iter/sec",
            "range": "stddev: 0.17112011880689199",
            "extra": "mean: 46.16666885359998 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_gauss]",
            "value": 0.6418551036283525,
            "unit": "iter/sec",
            "range": "stddev: 0.026055650511056752",
            "extra": "mean: 1.5579840284 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_uniform]",
            "value": 0.6559149505140279,
            "unit": "iter/sec",
            "range": "stddev: 0.014473674288897587",
            "extra": "mean: 1.5245879046000084 sec\nrounds: 5"
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
          "id": "8cb4e68f47dffa95209d8011c6ea242c07d7639f",
          "message": "1.1 ref sim ci workflow (#514)\n\n* added catalog file for 1.1 reference simulations to /data\r\n\r\n* possible approach to packaging reference simulation data files\r\n\r\n* moved files to better location\r\n\r\n* fixed up filepaths for 1.1_uniform reference simulation in /data\r\n\r\n* attempting to add benchmarks and updated ref_sim_1.1 in /data\r\n\r\n* benchmark needs more work but almost there\r\n\r\n* re-commented line\r\n\r\n* preliminary commit for benchmark workflow. Still need to set up 1.1_uniform for remote downloading.\r\n\r\n* pytest should now hopefully work remotely -- file on gdrive\r\n\r\n* trying slight change to see if all other tests run to completion properly\r\n\r\n* updated testsuite to install requests\r\n\r\n* test change\r\n\r\n* change and mpiexec\r\n\r\n* tried possibly fixing caching\r\n\r\n* benchmark action stuff\r\n\r\n* minor updates to see if alert works and if -np 4 speeds up run\r\n\r\n* test 3\r\n\r\n* just alert\r\n\r\n* test failed assert\r\n\r\n* swap back to ==\r\n\r\n* swapping back\r\n\r\n* small TODO\r\n\r\n* started adding new ref sim tests\r\n\r\n* formatting\r\n\r\n* added 1.1_gauss\r\n\r\n* got 1.1 uniform and gauss working, and resolved warnings\r\n\r\n* cosmetic update to testsuite\r\n\r\n* slight changes to test_run_ref.py\r\n\r\n* swapped to no longer writing the UVData objects and just returning one\r\n\r\n* changed gaussian beam to have proper syntax as well\r\n\r\n* preliminary attempt at sequential execution using matrix-lock\r\n\r\n* hopefully fixes the issue\r\n\r\n* 2nd attempt\r\n\r\n* had bad version of matrix-lock\r\n\r\n* removed matrix-lock\r\n\r\n* attempt to fix up caching\r\n\r\n* dummy commit\r\n\r\n* swapped back to triplicate for tests -- no longer doing command line input to pytest\r\n\r\n* added 1.3 sim files to data\r\n\r\n* added 1.2 simulations to data\r\n\r\n* google drive links\r\n\r\n* swapped workflow simulation run order\r\n\r\n* Swapped to downloading files from the BDR -- will trim workflow runtime down to ~1 hour using pedantic. Added mwa uvbeam sim files to data not yet tested.\r\n\r\n* figured out approach to parametrize all reference simulation tests using workflow. Still need to integrate mwa sims. Added pedantic benchmarking. Need to determine best approach to setting up workflow matrix.\r\n\r\n* filled out workflow refsim names, and added 1.1_mwa reference simulation to pytest and workflow.\r\n\r\n* changed the Brown Digital Repository file downloading to use a collection approach, added some print statements, added a line for formatting\r\n\r\n* removed 1.2_mwa files, minor comments change\r\n\r\n* Intermediate commit while attempting to switch benchmarking approach to use artifacts. Need to re-integrate the benchmark action, and create structure for concatenating benchmark output and uploading it.\r\n\r\n* fixed syntax error\r\n\r\n* commented out line to be re-added later\r\n\r\n* failed to comment out another line that skipped second part of workflow\r\n\r\n* test for python script which concatenates benchmarks\r\n\r\n* intermediate commit\r\n\r\n* first attempt at gh-pages approach\r\n\r\n* dummy change\r\n\r\n* preliminary approach to only pushing results of benchmark if push to 'main', and running only if Tests finishes successfully\r\n\r\n* removed dependence on Tests as that workflow seems to be failing independently\r\n\r\n* hopefully fixed yaml syntax\r\n\r\n* added initial output statistics to the reference simulation comparisons. currently only asserts '==' can implement others or even an absolute check\r\n\r\n* re-added setting history to be equal\r\n\r\n* fix\r\n\r\n* all current ref sims should run now, and implemented hopefully more robust downloading\r\n\r\n* commented out the 0 tolerance sim comparison check\r\n\r\n* added dummy counter (#513)\r\n\r\n* added dummy counter\r\n\r\n* [pre-commit.ci] auto fixes from pre-commit.com hooks\r\n\r\nfor more information, see https://pre-commit.ci\r\n\r\n---------\r\n\r\nCo-authored-by: pre-commit-ci[bot] <66853113+pre-commit-ci[bot]@users.noreply.github.com>\r\n\r\n* only one TODO left to resolve in test_run_ref (determining how strict object comparison should be)\r\n\r\n* cleaned up compare_post_benchmark.yaml a bit. now need to test running compare_post_benchmark using completion of another workflow (pull request and push)\r\n\r\n* updated approach to computing num_mismatched and fixed style\r\n\r\n* swapped compare_post_benchmark to run after Tests\r\n\r\n* minor edits to compare_post_benchmark so hopefully it runs\r\n\r\n* not sure why linking to tests isn't working -- swapping back\r\n\r\n* edited README / environment.yaml to discuss/require requests and pytest-benchmark\r\n\r\n* edited README\r\n\r\n* swapping to have defaults expected for pull request\r\n\r\n* changed underscore to hyphen to match style\r\n\r\n* Tentative README update -- should probably add a section of regression testing / ci  in the developer section of the docs, and amend README to link to it\r\n\r\n* made data comparison same as np.testing.assert_allclose defaults, removed some commented out code and comments\r\n\r\n* fixed typos in ci workflow\r\n\r\n* fixed formatting for a line\r\n\r\n* Futher updated the README\r\n\r\n* switching back to multiple ids\r\n\r\n* refactored job matrix\r\n\r\n* swapped discussion to docs for pytest regression testing\r\n\r\n---------\r\n\r\nCo-authored-by: pre-commit-ci[bot] <66853113+pre-commit-ci[bot]@users.noreply.github.com>",
          "timestamp": "2024-12-19T08:07:21-08:00",
          "tree_id": "f581b3ea4c0cb95a364a0bb8842db789a6ffd8c1",
          "url": "https://github.com/RadioAstronomySoftwareGroup/pyuvsim/commit/8cb4e68f47dffa95209d8011c6ea242c07d7639f"
        },
        "date": 1734625244587,
        "tool": "pytest",
        "benches": [
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_gauss]",
            "value": 0.0015316148227591836,
            "unit": "iter/sec",
            "range": "stddev: 0",
            "extra": "mean: 652.9056686710001 sec\nrounds: 1"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_uniform]",
            "value": 0.02199833779669115,
            "unit": "iter/sec",
            "range": "stddev: 0.4068254557319895",
            "extra": "mean: 45.45798001840001 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_uniform]",
            "value": 0.0015591259077841462,
            "unit": "iter/sec",
            "range": "stddev: 0",
            "extra": "mean: 641.3850190080001 sec\nrounds: 1"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_gauss]",
            "value": 0.021299215624968264,
            "unit": "iter/sec",
            "range": "stddev: 0.34236769288015534",
            "extra": "mean: 46.95008575000001 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_uniform]",
            "value": 0.6362727591334831,
            "unit": "iter/sec",
            "range": "stddev: 0.020287485070246032",
            "extra": "mean: 1.5716530146000025 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_mwa]",
            "value": 0.07819577129785138,
            "unit": "iter/sec",
            "range": "stddev: 0.11114088323826432",
            "extra": "mean: 12.788415324800019 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_gauss]",
            "value": 0.6447595720808491,
            "unit": "iter/sec",
            "range": "stddev: 0.01805891351413621",
            "extra": "mean: 1.5509657293999908 sec\nrounds: 5"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "66853113+pre-commit-ci[bot]@users.noreply.github.com",
            "name": "pre-commit-ci[bot]",
            "username": "pre-commit-ci[bot]"
          },
          "committer": {
            "email": "brynah@phys.washington.edu",
            "name": "Bryna Hazelton",
            "username": "bhazelton"
          },
          "distinct": true,
          "id": "d8253e99306b928fad9d231be272d5e2c36ee4e4",
          "message": "[pre-commit.ci] pre-commit autoupdate\n\nupdates:\n- [github.com/astral-sh/ruff-pre-commit: v0.8.3 → v0.8.6](https://github.com/astral-sh/ruff-pre-commit/compare/v0.8.3...v0.8.6)",
          "timestamp": "2025-01-07T10:36:57-07:00",
          "tree_id": "7ae38b9a30b3805247158d95c3896ff9d92e4ce9",
          "url": "https://github.com/RadioAstronomySoftwareGroup/pyuvsim/commit/d8253e99306b928fad9d231be272d5e2c36ee4e4"
        },
        "date": 1736272205130,
        "tool": "pytest",
        "benches": [
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_gauss]",
            "value": 0.0015521162975521727,
            "unit": "iter/sec",
            "range": "stddev: 0",
            "extra": "mean: 644.2816183150001 sec\nrounds: 1"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_uniform]",
            "value": 0.021727920687523387,
            "unit": "iter/sec",
            "range": "stddev: 0.2716723886891983",
            "extra": "mean: 46.023732062600004 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_uniform]",
            "value": 0.001526862590887347,
            "unit": "iter/sec",
            "range": "stddev: 0",
            "extra": "mean: 654.937782855 sec\nrounds: 1"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_gauss]",
            "value": 0.021597186964759908,
            "unit": "iter/sec",
            "range": "stddev: 0.07158112960252856",
            "extra": "mean: 46.302326392400005 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_uniform]",
            "value": 0.6312704993526635,
            "unit": "iter/sec",
            "range": "stddev: 0.03003446692160079",
            "extra": "mean: 1.5841069731999995 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_mwa]",
            "value": 0.07892553528927293,
            "unit": "iter/sec",
            "range": "stddev: 0.09678511045031156",
            "extra": "mean: 12.670170640399999 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_gauss]",
            "value": 0.6685658909425418,
            "unit": "iter/sec",
            "range": "stddev: 0.021678772385816808",
            "extra": "mean: 1.4957388845999957 sec\nrounds: 5"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "bryna.hazelton@gmail.com",
            "name": "Bryna Hazelton",
            "username": "bhazelton"
          },
          "committer": {
            "email": "brynah@phys.washington.edu",
            "name": "Bryna Hazelton",
            "username": "bhazelton"
          },
          "distinct": true,
          "id": "1d89b554ee34c7ac8dbda177cb2c841fe597048b",
          "message": "fix missing git tag info",
          "timestamp": "2025-01-16T12:18:37-08:00",
          "tree_id": "60b235ff309df2ae40353e75a6acd4b679d9eef2",
          "url": "https://github.com/RadioAstronomySoftwareGroup/pyuvsim/commit/1d89b554ee34c7ac8dbda177cb2c841fe597048b"
        },
        "date": 1737059528070,
        "tool": "pytest",
        "benches": [
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_gauss]",
            "value": 0.001560339523469184,
            "unit": "iter/sec",
            "range": "stddev: 0",
            "extra": "mean: 640.88615648 sec\nrounds: 1"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_uniform]",
            "value": 0.021596131890939667,
            "unit": "iter/sec",
            "range": "stddev: 0.5287187286860251",
            "extra": "mean: 46.30458848139999 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_uniform]",
            "value": 0.0015699655582524538,
            "unit": "iter/sec",
            "range": "stddev: 0",
            "extra": "mean: 636.9566483440001 sec\nrounds: 1"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_gauss]",
            "value": 0.02103319212483775,
            "unit": "iter/sec",
            "range": "stddev: 0.27272586895706724",
            "extra": "mean: 47.54390080520001 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_uniform]",
            "value": 0.6262136818010258,
            "unit": "iter/sec",
            "range": "stddev: 0.0326080377861726",
            "extra": "mean: 1.596898996399989 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_mwa]",
            "value": 0.07868020119813834,
            "unit": "iter/sec",
            "range": "stddev: 0.03731063591471391",
            "extra": "mean: 12.709677717800002 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_gauss]",
            "value": 0.6271819227906298,
            "unit": "iter/sec",
            "range": "stddev: 0.03151050295941445",
            "extra": "mean: 1.5944337100000041 sec\nrounds: 5"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "66853113+pre-commit-ci[bot]@users.noreply.github.com",
            "name": "pre-commit-ci[bot]",
            "username": "pre-commit-ci[bot]"
          },
          "committer": {
            "email": "brynah@phys.washington.edu",
            "name": "Bryna Hazelton",
            "username": "bhazelton"
          },
          "distinct": true,
          "id": "118c21b66b87bf82ca94f3c9fff2a0a7c338269e",
          "message": "[pre-commit.ci] auto fixes from pre-commit.com hooks\n\nfor more information, see https://pre-commit.ci",
          "timestamp": "2025-01-17T09:39:36-08:00",
          "tree_id": "af8befc46a05d9374e2248a947efa8b8b5c30849",
          "url": "https://github.com/RadioAstronomySoftwareGroup/pyuvsim/commit/118c21b66b87bf82ca94f3c9fff2a0a7c338269e"
        },
        "date": 1737136359796,
        "tool": "pytest",
        "benches": [
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_gauss]",
            "value": 0.0015479437756907336,
            "unit": "iter/sec",
            "range": "stddev: 0",
            "extra": "mean: 646.018295822 sec\nrounds: 1"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_uniform]",
            "value": 0.02177942538886528,
            "unit": "iter/sec",
            "range": "stddev: 0.19159792490829375",
            "extra": "mean: 45.914893627599994 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_uniform]",
            "value": 0.0015603959932725341,
            "unit": "iter/sec",
            "range": "stddev: 0",
            "extra": "mean: 640.8629631910001 sec\nrounds: 1"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_gauss]",
            "value": 0.02135584114686689,
            "unit": "iter/sec",
            "range": "stddev: 0.19160984977847983",
            "extra": "mean: 46.825596478399994 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_uniform]",
            "value": 0.6370955056870747,
            "unit": "iter/sec",
            "range": "stddev: 0.01944800202413224",
            "extra": "mean: 1.5696233783999958 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_mwa]",
            "value": 0.07903335325480444,
            "unit": "iter/sec",
            "range": "stddev: 0.05781118816419168",
            "extra": "mean: 12.652885886999991 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_gauss]",
            "value": 0.6518087670799065,
            "unit": "iter/sec",
            "range": "stddev: 0.017021689862788392",
            "extra": "mean: 1.5341923129999997 sec\nrounds: 5"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "66853113+pre-commit-ci[bot]@users.noreply.github.com",
            "name": "pre-commit-ci[bot]",
            "username": "pre-commit-ci[bot]"
          },
          "committer": {
            "email": "brynah@phys.washington.edu",
            "name": "Bryna Hazelton",
            "username": "bhazelton"
          },
          "distinct": true,
          "id": "c93c5078e92b49fcd823c136ecdc3dcf7e198517",
          "message": "[pre-commit.ci] pre-commit autoupdate\n\nupdates:\n- [github.com/astral-sh/ruff-pre-commit: v0.9.1 → v0.9.2](https://github.com/astral-sh/ruff-pre-commit/compare/v0.9.1...v0.9.2)",
          "timestamp": "2025-01-21T10:01:09-08:00",
          "tree_id": "1fdcbd65dd393f6387eca7d5137a6e006d48094f",
          "url": "https://github.com/RadioAstronomySoftwareGroup/pyuvsim/commit/c93c5078e92b49fcd823c136ecdc3dcf7e198517"
        },
        "date": 1737483254203,
        "tool": "pytest",
        "benches": [
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_gauss]",
            "value": 0.0015769136267966106,
            "unit": "iter/sec",
            "range": "stddev: 0",
            "extra": "mean: 634.150141775 sec\nrounds: 1"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_uniform]",
            "value": 0.021694298808433373,
            "unit": "iter/sec",
            "range": "stddev: 0.09426544722809276",
            "extra": "mean: 46.0950597588 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_uniform]",
            "value": 0.0015591317844578009,
            "unit": "iter/sec",
            "range": "stddev: 0",
            "extra": "mean: 641.382601502 sec\nrounds: 1"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_gauss]",
            "value": 0.022066760935207385,
            "unit": "iter/sec",
            "range": "stddev: 0.22351552058958718",
            "extra": "mean: 45.31702694999999 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_uniform]",
            "value": 0.642652685117938,
            "unit": "iter/sec",
            "range": "stddev: 0.023043254318960975",
            "extra": "mean: 1.55605045019999 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_mwa]",
            "value": 0.07641653919470928,
            "unit": "iter/sec",
            "range": "stddev: 0.12416254309554707",
            "extra": "mean: 13.086172320000003 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_gauss]",
            "value": 0.6301092142567387,
            "unit": "iter/sec",
            "range": "stddev: 0.023950054673965035",
            "extra": "mean: 1.587026466800006 sec\nrounds: 5"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "66853113+pre-commit-ci[bot]@users.noreply.github.com",
            "name": "pre-commit-ci[bot]",
            "username": "pre-commit-ci[bot]"
          },
          "committer": {
            "email": "brynah@phys.washington.edu",
            "name": "Bryna Hazelton",
            "username": "bhazelton"
          },
          "distinct": true,
          "id": "2b36b9f11e3d70ac937e24b97d5078d256683bba",
          "message": "[pre-commit.ci] pre-commit autoupdate\n\nupdates:\n- [github.com/astral-sh/ruff-pre-commit: v0.9.2 → v0.9.3](https://github.com/astral-sh/ruff-pre-commit/compare/v0.9.2...v0.9.3)",
          "timestamp": "2025-01-28T10:34:41-08:00",
          "tree_id": "5fa72b84134468857561fcb8028ee8572c327466",
          "url": "https://github.com/RadioAstronomySoftwareGroup/pyuvsim/commit/2b36b9f11e3d70ac937e24b97d5078d256683bba"
        },
        "date": 1738090167314,
        "tool": "pytest",
        "benches": [
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_gauss]",
            "value": 0.6454457602559188,
            "unit": "iter/sec",
            "range": "stddev: 0.01924205514411844",
            "extra": "mean: 1.5493168622000097 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_uniform]",
            "value": 0.6384025192270474,
            "unit": "iter/sec",
            "range": "stddev: 0.02092669822326019",
            "extra": "mean: 1.5664098587999944 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_uniform]",
            "value": 0.0015646581864819157,
            "unit": "iter/sec",
            "range": "stddev: 0",
            "extra": "mean: 639.117226139 sec\nrounds: 1"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_gauss]",
            "value": 0.021950290937387876,
            "unit": "iter/sec",
            "range": "stddev: 0.14509486896487162",
            "extra": "mean: 45.55748271640001 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.2_gauss]",
            "value": 0.0015440168248602503,
            "unit": "iter/sec",
            "range": "stddev: 0",
            "extra": "mean: 647.661336262 sec\nrounds: 1"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.3_uniform]",
            "value": 0.02138753171448189,
            "unit": "iter/sec",
            "range": "stddev: 0.5134439862560964",
            "extra": "mean: 46.75621354299999 sec\nrounds: 5"
          },
          {
            "name": "tests/test_run_ref.py::test_run_sim[1.1_mwa]",
            "value": 0.07917425603763331,
            "unit": "iter/sec",
            "range": "stddev: 0.09452392079117607",
            "extra": "mean: 12.630368127799994 sec\nrounds: 5"
          }
        ]
      }
    ]
  }
}