"""datalayer/autocache.py: the @cacheFunction decorator and its helpers.

Project is patched at the `hera` package level (not autocache's own
namespace, since it imports Project lazily inside each function) to route
through unit_files_directory instead of the guarded fake home.
"""
import pytest

from hera.datalayer import autocache


@pytest.fixture(autouse=True)
def _patched_project(monkeypatch, unit_files_directory):
    import hera
    from hera.datalayer import Project as RealProject

    def _patched(projectName=None, **kwargs):
        kwargs.setdefault("filesDirectory", unit_files_directory)
        return RealProject(projectName=projectName, **kwargs)

    monkeypatch.setattr(hera, "Project", _patched)


@pytest.mark.unit
class TestSerializationHelpers:
    def test_is_mongo_serializable_true_for_plain_data(self):
        assert autocache.cacheDecorators.is_mongo_serializable(1) is True

    def test_is_mongo_serializable_false_for_a_function(self):
        assert autocache.cacheDecorators.is_mongo_serializable(lambda x: x) is False

    def test_obj_to_txt_round_trips(self):
        txt = autocache.cacheDecorators.obj_to_txt({"a": 1})
        assert autocache.cacheDecorators.txt_to_obj(txt) == {"a": 1}


@pytest.mark.unit
class TestGetFullFuncName:
    def test_plain_function_uses_qualname(self):
        def foo(x):
            pass

        cd = autocache.cacheDecorators(func=foo, dataFormat=None)
        assert cd._get_full_func_name(foo) == foo.__qualname__

    def test_bound_method_includes_class_and_module(self):
        class C:
            def bar(self):
                pass

        cd = autocache.cacheDecorators(func=C().bar, dataFormat=None)
        name = cd._get_full_func_name(C().bar)
        assert name.endswith("C.bar")

    def test_non_callable_raises(self):
        cd = autocache.cacheDecorators(func=lambda: None, dataFormat=None)
        with pytest.raises(TypeError, match="not callable"):
            cd._get_full_func_name(42)


@pytest.mark.unit
class TestCacheFunctionDecorator:
    def test_second_call_uses_the_cache_not_a_recompute(self):
        calls = []

        @autocache.cacheFunction(projectName="CACHE_TEST_1")
        def compute(x):
            calls.append(x)
            return x * 2

        assert compute(5) == 10
        assert compute(5) == 10
        assert calls == [5]

    def test_different_arguments_are_cached_separately(self):
        calls = []

        @autocache.cacheFunction(projectName="CACHE_TEST_2")
        def compute(x):
            calls.append(x)
            return x * 2

        compute(1)
        compute(2)
        assert calls == [1, 2]

    def test_clear_function_cache_forces_a_recompute(self):
        calls = []

        @autocache.cacheFunction(projectName="CACHE_TEST_3")
        def compute(x):
            calls.append(x)
            return x * 2

        compute(5)
        autocache.clearFunctionCache(functionName=None, projectName="CACHE_TEST_3")
        compute(5)
        assert calls == [5, 5]

    def test_clear_all_functions_cache_delegates_to_clear_function_cache(self):
        calls = []

        @autocache.cacheFunction(projectName="CACHE_TEST_4")
        def compute(x):
            calls.append(x)
            return x * 2

        compute(5)
        autocache.clearAllFunctionsCache(projectName="CACHE_TEST_4")
        compute(5)
        assert calls == [5, 5]

    def test_return_doc_true_returns_a_tuple(self):
        @autocache.cacheFunction(projectName="CACHE_TEST_5", returnDoc=True)
        def compute(x):
            return x * 2

        result, doc = compute(5)
        assert result == 10
        assert doc is not None

    def test_post_process_function_transforms_the_result(self):
        @autocache.cacheFunction(projectName="CACHE_TEST_6", postProcessFunction=lambda r: r + 1)
        def compute(x):
            return x * 2

        assert compute(5) == 11
