"""respt 测试包。

运行全部测试::

    python -m unittest discover -s tests -t . -v

只跑单元测试（不需要外部工具与数据库）::

    python -m unittest tests.test_intervals tests.test_blast tests.test_annotate \
        tests.test_context tests.test_utils tests.test_cli -v

端到端测试（真实数据库 + stub 外部工具）::

    python -m unittest tests.test_e2e -v
"""
