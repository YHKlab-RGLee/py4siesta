"""One provider-independent chat-model initialization and validation boundary."""

import importlib.util
import json
import os
from typing import Any, Dict

from pydantic import ValidationError

from .types import AgentConfigurationError, ModelOutputError


SUPPORTED_PROVIDERS = {
    "openai": ("langchain_openai", "OPENAI_API_KEY"),
    "anthropic": ("langchain_anthropic", "ANTHROPIC_API_KEY"),
    "ollama": ("langchain_ollama", None),
}


class ChatModelClient:
    """Hide LangChain/provider response objects from the rest of the package."""

    def __init__(self, model):
        self._model = model

    def structured(self, instructions, user_input, schema):
        prompt = [
            ("system", instructions),
            ("human", user_input),
        ]
        try:
            result = self._model.with_structured_output(schema).invoke(prompt)
        except (AttributeError, NotImplementedError):
            result = self._json_fallback(instructions, user_input, schema)
        except Exception as exc:
            raise AgentConfigurationError(
                "The configured model could not produce structured output: %s" % exc
            ) from exc

        try:
            if isinstance(result, schema):
                return result
            if isinstance(result, dict):
                return schema.model_validate(result)
            return schema.model_validate_json(result)
        except (TypeError, ValidationError, ValueError) as exc:
            raise ModelOutputError("Model structured output failed validation: %s" % exc) from exc

    def _json_fallback(self, instructions, user_input, schema):
        json_schema = json.dumps(schema.model_json_schema(), sort_keys=True)
        response = self._model.invoke(
            [
                (
                    "system",
                    instructions
                    + "\nReturn only JSON matching this schema exactly:\n"
                    + json_schema,
                ),
                ("human", user_input),
            ]
        )
        content = getattr(response, "content", None)
        if not isinstance(content, str):
            raise AgentConfigurationError(
                "The configured model supports neither native structured output "
                "nor a plain-text JSON fallback."
            )
        return content


def create_model_client(
    provider=None,
    model=None,
    environment: Dict[str, str] = None,
    initializer=None,
):
    """Create the shared model client from ``PY4SIESTA_LLM_*`` configuration."""

    environment = os.environ if environment is None else environment
    provider = (provider or environment.get("PY4SIESTA_LLM_PROVIDER", "")).lower()
    model = model or environment.get("PY4SIESTA_LLM_MODEL")
    if not provider or not model:
        raise AgentConfigurationError(
            "Set PY4SIESTA_LLM_PROVIDER and PY4SIESTA_LLM_MODEL."
        )
    if provider not in SUPPORTED_PROVIDERS:
        raise AgentConfigurationError(
            "Unsupported model provider %r. Supported providers: %s."
            % (provider, ", ".join(sorted(SUPPORTED_PROVIDERS)))
        )

    package, credential = SUPPORTED_PROVIDERS[provider]
    if importlib.util.find_spec(package) is None:
        raise AgentConfigurationError(
            "Provider %s requires the optional package %s." % (provider, package)
        )
    if credential and not environment.get(credential):
        raise AgentConfigurationError(
            "%s is required for provider %s." % (credential, provider)
        )

    if initializer is None:
        from langchain.chat_models import init_chat_model

        initializer = init_chat_model
    try:
        chat_model = initializer(model=model, model_provider=provider, temperature=0)
    except Exception as exc:
        raise AgentConfigurationError(
            "Could not initialize provider %s model %s: %s" % (provider, model, exc)
        ) from exc
    return ChatModelClient(chat_model)
