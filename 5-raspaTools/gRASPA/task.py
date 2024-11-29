import BasicSet, Framework, ConditionSet, MoviesSet, Component
from BasicSet import DEFAULT_SETTINGS as BASIC_DEFAULTS
from Framework import DEFAULT_SETTINGS as FRAMEWORK_DEFAULTS
from ConditionSet import DEFAULT_SETTINGS as CONDITIION_DEFAULTS
from MoviesSet import DEFAULT_SETTINGS as MOVIES_DEFAULTS
from Component import DEFAULT_SETTINGS as COMPONENT_DEFAULTS


class TaskManager:
    _default_settings = {}

    @classmethod
    def register_defaults(cls, task_name, defaults):
        if callable(defaults):
            cls._default_settings[task_name] = cls.extract_defaults_from_function(defaults)
        elif isinstance(defaults, dict):
            cls._default_settings[task_name] = defaults.copy()
        else:
            raise TypeError(f"Defaults for task '{task_name}' must be a callable or a dictionary.")

    @classmethod
    def extract_defaults_from_function(cls, func):
        import inspect
        signature = inspect.signature(func)
        return {
            key: value.default
            for key, value in signature.parameters.items()
            if value.default is not inspect.Parameter.empty
        }

    @classmethod
    def get_defaults(cls, task_name):
        if task_name not in cls._default_settings:
            raise ValueError(f"check your '{task_name}'.")
        return cls._default_settings[task_name]

    @classmethod
    def update_defaults(cls, task_name, **kwargs):
        if task_name not in cls._default_settings:
            raise ValueError(f"check your '{task_name}'.")
        cls._default_settings[task_name].update(kwargs)


class Generate:
    MODULE_MAPPING = {
                        "BasicSet": BasicSet.Generate,
                        "Framework": Framework.Generate,
                        "ConditionSet": ConditionSet.Generate,
                        "MoviesSet": MoviesSet.Generate,
                        "Component": Component.Generate,
                    }

    def __init__(self, task, out_folder, **kwargs):
        print(f"user settings: {kwargs}")
        if task not in ["SingleAdsorption", "MixtureAdsorption"]:
            raise ValueError("Invalid task name. Please check your task name.")
        self.out_folder = out_folder
        self.task = task
        self.kwargs = kwargs
        if "Component" in self.kwargs:
            self.kwargs["Component"]["task"] = self.task
        self.result = self._generate_config()
        with open(f"{self.out_folder}/simulation.input", "w") as file:
            file.write(self.result)

    def _generate_config(self):
        structure = self.kwargs.get("structure")
        if not structure:
            raise ValueError("Structure file is not specified.")
        task_kwargs = self.kwargs.get("task", {})
        task_kwargs["task"] = self.task
        framework_kwargs = self.kwargs.get("Framework", {})
        framework_kwargs["structure"] = structure
        config_parts = []
        for module_name, generate_func in self.MODULE_MAPPING.items():
            module_kwargs = self.kwargs.get(module_name, {})
            if module_name == "Framework":
                module_kwargs = framework_kwargs
            if module_name == "task":
                module_kwargs = task_kwargs
            instance = generate_func(**module_kwargs)
            config_parts.append(instance.generate_config())
        return "\n".join(config_parts)
    
